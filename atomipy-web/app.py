import os
import sys
import uuid
import zipfile
import io
import time
import logging
import traceback
import importlib.util
import copy
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash, session, Response, jsonify
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge
from flask_executor import Executor
import shutil

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('atomipy-web')

# Global settings for timeouts and limits
PROCESSING_TIMEOUT = 600  # seconds
MAX_FILE_SIZE = 32 * 1024 * 1024  # 32 MB
DEFAULT_CIF_FUSE_RMAX = 0.85
DEFAULT_CIF_FUSE_CRITERIA = "average"
DEFAULT_BVS_DELTA_THRESHOLD = -0.50
DEFAULT_BVS_MAX_ADDITIONS = 10
ANGLE_TERM_OPTIONS = {"none", "0", "250", "500", "1500"}
DEFAULT_ANGLE_TERM_MINFF = "500"
DEFAULT_ANGLE_TERM_CLAYFF = "none"

# Add the parent directory to sys.path to import atomipy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Lazy loader for atomipy to reduce initial memory footprint
_ap = None
def get_ap():
    global _ap
    if _ap is None:
        import atomipy
        _ap = atomipy
    return _ap

# Direct lazy getters for frequently used submodules if needed
def get_ap_forcefield():
    ap = get_ap()
    from atomipy import forcefield
    return forcefield

# Base directory for resolving templates/static regardless of CWD
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

app = Flask(
    __name__,
    template_folder=os.path.join(BASE_DIR, 'templates'),
    static_folder=os.path.join(BASE_DIR, 'static')
)
app.config['SECRET_KEY'] = os.urandom(24)
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE  # Reduced to 16 MB for better stability

# App Engine specific configuration
app.config['EXECUTOR_TYPE'] = 'thread'
app.config['EXECUTOR_MAX_WORKERS'] = 8  # Limit concurrent tasks to prevent resource exhaustion

# Initialize Flask-Executor
executor = Executor(app)
tasks_status = {}  # Dictionary to store task status

# Define upload and results folders
UPLOAD_FOLDER_NAME = 'uploads'
RESULTS_FOLDER_NAME = 'results'

# Use /tmp for writable storage on App Engine, otherwise use local folders
# Check for any Google Cloud environment: GAE_ENV, GAE_INSTANCE, or CLOUD_RUN_JOB
if os.environ.get('GAE_ENV', '').startswith('standard') or os.environ.get('GAE_INSTANCE') or os.environ.get('CLOUD_RUN_JOB'):
    app.config['UPLOAD_FOLDER'] = '/tmp/uploads'
    app.config['RESULTS_FOLDER'] = '/tmp/results'
    app.config['TEMP_DIR'] = '/tmp'
    print(f"Running in cloud environment, using /tmp for writable storage")
else:
    app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), UPLOAD_FOLDER_NAME)
    app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), RESULTS_FOLDER_NAME)

app.config['ALLOWED_EXTENSIONS'] = {'gro', 'pdb', 'xyz', 'cif', 'mmcif', 'mcif'}


def gemmi_is_available():
    return importlib.util.find_spec("gemmi") is not None


app.config['GEMMI_AVAILABLE'] = gemmi_is_available()


def parse_bool(value, default=False):
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def parse_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def add_hydrogens_from_bvs(
    atoms,
    box_dim,
    delta_threshold=DEFAULT_BVS_DELTA_THRESHOLD,
    max_additions=DEFAULT_BVS_MAX_ADDITIONS,
):
    """
    Add one H to strongly underbonded oxygen sites selected via BVS.

    Returns
    -------
    tuple
        (updated_atoms, added_h_count, selected_site_count)
    """
    if max_additions <= 0:
        return atoms, 0, 0

    bvs_report = get_ap().analyze_bvs(copy.deepcopy(atoms), box_dim, top_n=max_additions)
    bvs_results = bvs_report.get("results", [])

    candidates = []
    for row in bvs_results:
        element = (row.get("element") or "").upper()
        delta = row.get("delta")
        if element != "O" or delta is None or delta >= delta_threshold:
            continue

        idx0 = int(row.get("index", 0)) - 1
        if idx0 < 0 or idx0 >= len(atoms):
            continue

        # Skip oxygens that are already bonded to at least one hydrogen.
        bonded_h = False
        for neigh_idx, *_ in row.get("bonds", []):
            j = int(neigh_idx) - 1
            if 0 <= j < len(atoms):
                neigh_el = (atoms[j].get("element") or "").upper()
                if neigh_el == "H":
                    bonded_h = True
                    break
        if bonded_h:
            continue

        candidates.append((delta, idx0, row))

    candidates.sort(key=lambda item: item[0])  # most negative delta first
    candidates = candidates[:max_additions]

    added_h = 0
    for _, idx0, row in candidates:
        if idx0 >= len(atoms):
            continue

        original_type = atoms[idx0].get("type", "O")
        marker_type = f"__BVS_TARGET_{idx0}__"
        atoms[idx0]["type"] = marker_type

        prev_count = len(atoms)
        cn = row.get("cn")
        coordination_target = max(1, int(cn) + 1) if isinstance(cn, (int, float)) else 2
        ap = get_ap()
        atoms = ap.add_H_atom(
            atoms,
            box_dim,
            target_type=marker_type,
            h_type="H",
            bond_length=0.98,
            coordination=coordination_target,
            max_h_per_atom=1,
        )

        if idx0 < len(atoms):
            atoms[idx0]["type"] = original_type
        added_h += max(0, len(atoms) - prev_count)

    return atoms, added_h, len(candidates)

# Create upload and results folders if they don't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

# --- Background Task ---
def process_file_task(
    task_id,
    filepath,
    filename,
    ff_type,
    output_formats,
    results_id,
    results_dir,
    generate_topology=True,
    expand_symmetry=True,
    fuse_overlaps=True,
    fuse_rmax=DEFAULT_CIF_FUSE_RMAX,
    fuse_criteria=DEFAULT_CIF_FUSE_CRITERIA,
    add_hydrogen_bvs=False,
    bvs_delta_threshold=DEFAULT_BVS_DELTA_THRESHOLD,
    bvs_max_additions=DEFAULT_BVS_MAX_ADDITIONS,
    replicate_structure=False,
    replicate_nx=1,
    replicate_ny=1,
    replicate_nz=1,
    angle_terms=DEFAULT_ANGLE_TERM_MINFF,
    reset_molid=True,
):
    start_time = time.time()
    logger.info(f"Starting task {task_id} for file {filename}")
    base_filename = filename.rsplit('.', 1)[0]
    try:
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Reading structure', 'progress': 10}
        # Use appropriate import function based on file extension
        file_extension = filename.rsplit('.', 1)[1].lower()
        atoms, Cell, Box_dim = None, None, None
        ap = get_ap()
        if file_extension == 'gro':
            atoms, Box_dim = ap.import_gro(filepath)
            if Box_dim is None:
                raise ValueError("GRO file is missing box dimensions in the last line.")
            Cell = ap.Box_dim2Cell(Box_dim)
        elif file_extension == 'pdb':
            atoms, Cell = ap.import_pdb(filepath)
            if Cell is None:
                raise ValueError("PDB file must contain box dimensions in a CRYST1 record.")
            Box_dim = ap.Cell2Box_dim(Cell)
        elif file_extension == 'xyz':
            # Import XYZ file with Box dimensions from the comment line
            atoms, Cell = ap.import_xyz(filepath)
            if Cell is None:
                raise ValueError("XYZ file must contain Box dimensions on the second line after a # character (e.g., # 40.0 40.0 40.0)")
            Box_dim = ap.Cell2Box_dim(Cell)
        elif file_extension in {'cif', 'mmcif', 'mcif'}:
            atoms, Cell = ap.import_cif(filepath, expand_symmetry=expand_symmetry)
            if Cell is None:
                raise ValueError("CIF file is missing unit cell parameters.")
            Box_dim = ap.Cell2Box_dim(Cell)
        else:
            raise ValueError(f"Unsupported file extension: {file_extension}")

        if fuse_overlaps:
            tasks_status[task_id] = {
                'status': 'Processing',
                'step': f'Fusing overlapping sites (r={fuse_rmax:.2f} Å)',
                'progress': 18,
            }
            n_before_fuse = len(atoms)
            atoms = ap.fuse_atoms(atoms, Box_dim, rmax=fuse_rmax, criteria=fuse_criteria)
            logger.info(
                "Fuse step (%s): %d -> %d atoms (criteria=%s, rmax=%.3f)",
                filename,
                n_before_fuse,
                len(atoms),
                fuse_criteria,
                fuse_rmax,
            )

        if add_hydrogen_bvs:
            tasks_status[task_id] = {
                'status': 'Processing',
                'step': 'Protonating underbonded oxygen sites (BVS)',
                'progress': 22,
            }
            atoms, added_h, selected_sites = add_hydrogens_from_bvs(
                atoms,
                Box_dim,
                delta_threshold=bvs_delta_threshold,
                max_additions=bvs_max_additions,
            )
            logger.info(
                "BVS H-addition step (%s): selected_sites=%d, added_H=%d, delta_threshold=%.3f",
                filename,
                selected_sites,
                added_h,
                bvs_delta_threshold,
            )

        replicate_factors = [
            max(1, parse_int(replicate_nx, 1)),
            max(1, parse_int(replicate_ny, 1)),
            max(1, parse_int(replicate_nz, 1)),
        ]
        if replicate_structure and replicate_factors != [1, 1, 1]:
            tasks_status[task_id] = {
                'status': 'Processing',
                'step': f'Replicating structure ({replicate_factors[0]}x{replicate_factors[1]}x{replicate_factors[2]})',
                'progress': 26,
            }
            n_before_replicate = len(atoms)
            atoms, Box_dim, Cell = ap.replicate_system(
                atoms,
                Box=Box_dim,
                replicate=replicate_factors,
            )
            logger.info(
                "Replicate step (%s): %d -> %d atoms (replicate=%s)",
                filename,
                n_before_replicate,
                len(atoms),
                replicate_factors,
            )

        tasks_status[task_id] = {'status': 'Processing', 'step': f'Assigning {ff_type} atom types', 'progress': 30}
        if ff_type == 'minff':
            # Import necessary functions from atomipy
            from atomipy.forcefield import minff

            # Separate water from mineral+ions BEFORE running MINFF so that water
            # oxygens do not perturb the coordination-number detection for Al/Mg.
            # This mirrors the run_sys2minff.py pipeline exactly.
            # Only do this when reset_molid is True (the default / recommended option).
            if reset_molid:
                try:
                    ap = get_ap()
                    SOL, noSOL = ap.find_H2O(atoms, Box_dim)
                    noSOL = ap.assign_resname(noSOL)
                    MIN = [a for a in noSOL if a.get('resname') == 'MIN']
                    OTHER = [a for a in noSOL if a.get('resname') != 'MIN']
                    if MIN:
                        MIN = ap.update(MIN, molid=1)
                    atoms = ap.update(OTHER, MIN, SOL)
                except Exception as e:
                    print(f"Warning: water/ion separation failed ({e}); running minff on full system.")

            # Generate log file path in the writable results directory
            log_file = os.path.join(results_dir, 'minff_structure_stats.log')

            # Assign atom types using MINFF forcefield and generate statistics in one call
            atoms = get_ap_forcefield().minff(atoms, Box=Box_dim, log=True, log_file=log_file)
        elif ff_type == 'clayff':
            # Generate log file path in the writable results directory
            log_file = os.path.join(results_dir, 'clayff_structure_stats.log')

            # Assign atom types using CLAYFF forcefield and generate statistics in one call
            atoms = get_ap_forcefield().clayff(atoms, Box=Box_dim, log=True, log_file=log_file)
        tasks_status[task_id] = {'status': 'Processing', 'step': f'Calculating charges ({ff_type})', 'progress': 50}
        # Comprehensive debug of the structure
        print(f"Type of atoms after processing: {type(atoms)}")
        
        # Special handling for tuple returned from forcefield functions
        if isinstance(atoms, tuple) and len(atoms) == 2:
            # If atoms is a tuple of (atoms_list, Box_dim), extract the atoms list
            atoms_from_tuple, Box_dim_from_tuple = atoms
            atoms = atoms_from_tuple
            print(f"Extracted atoms list from tuple, new Box_dim: {Box_dim_from_tuple}")
            # Update Box_dim if needed
            if Box_dim_from_tuple is not None and len(Box_dim_from_tuple) in [3, 9]:
                Box_dim = Box_dim_from_tuple
                print(f"Updated Box_dim to: {Box_dim}")
        
        # Ensure we have a valid list to work with
        if not isinstance(atoms, list):
            print(f"WARNING: atoms is not a list but {type(atoms)}")
            if atoms is None:
                atoms = []
            elif hasattr(atoms, '__iter__'):
                atoms = list(atoms)
            else:
                atoms = [atoms] if atoms else []
            print(f"Converted atoms to a list with {len(atoms)} items")
        if atoms and len(atoms) > 0:
            print(f"Type of first atom: {type(atoms[0])}")
            print(f"Sample atom contents: {atoms[0]}")
            # Create a new list for the converted atoms
            new_atoms = []
            conversion_count = 0
            for i, atom in enumerate(atoms):
                # For debugging, print every 1000th atom
                if i % 1000 == 0:
                    print(f"Processing atom {i}, type: {type(atom)}")
                if not isinstance(atom, dict):
                    conversion_count += 1
                    # Handle the case where atoms are lists or tuples
                    if isinstance(atom, (list, tuple)):
                        # Use common atom properties, adjust as needed based on atom fields
                        atom_dict = {}
                        # Try to infer property names from length of tuple
                        if len(atom) <= 12:  # Common case for molecular data
                            properties = ['molid', 'index', 'resname', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'type', 'charge', 'mass']
                            for j, prop in enumerate(properties):
                                if j < len(atom):
                                    atom_dict[prop] = atom[j]
                        else:
                            # Generic case for longer tuples
                            for j, value in enumerate(atom):
                                atom_dict[f'prop{j}'] = value
                        new_atoms.append(atom_dict)
                    else:
                        print(f"WARNING: Cannot convert atom {i} of type {type(atom)} to dictionary")
                        # Add a placeholder to maintain index consistency
                        new_atoms.append({'index': i, 'warning': f'Invalid type: {type(atom)}'})
                else:
                    new_atoms.append(atom)
            if conversion_count > 0:
                print(f"Converted {conversion_count} atoms to dictionary format")
                atoms = new_atoms
        # Configure angle treatment for topology generation.
        include_angles = angle_terms != "none"
        topology_kangle = 0 if angle_terms == "none" else parse_int(angle_terms, 500)
        topology_max_angle = None if include_angles else 0.0

        # Generate selected output files
        generated_files = []
        progress_step = 70
        
        # Only generate topology files if requested
        if generate_topology and output_formats:
            progress_increment = (85 - progress_step) / len(output_formats) if output_formats else 0
            if 'itp' in output_formats:
                tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing ITP', 'progress': int(progress_step)}
                topology_itp = os.path.join(results_dir, f"{base_filename}_{ff_type}.itp")
                from atomipy import write_top
                write_top.itp(
                    atoms,
                    Box=Box_dim,
                    file_path=topology_itp,
                    explicit_angles=1 if include_angles else 0,
                    KANGLE=topology_kangle,
                    max_angle=topology_max_angle,
                )
                generated_files.append(topology_itp)
                print(f"ITP file for {ff_type} written successfully to {topology_itp}")
                progress_step += progress_increment
            if 'psf' in output_formats:
                tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing PSF', 'progress': int(progress_step)}
                topology_psf = os.path.join(results_dir, f"{base_filename}_{ff_type}.psf")
                from atomipy import write_top
                write_top.psf(
                    atoms,
                    Box=Box_dim,
                    file_path=topology_psf,
                    max_angle=topology_max_angle,
                )
                generated_files.append(topology_psf)
                print(f"PSF file for {ff_type} written successfully to {topology_psf}")
                progress_step += progress_increment
            if 'lmp' in output_formats:  # Check for 'lmp' from form value
                tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing LAMMPS Data', 'progress': int(progress_step)}
                topology_lmp = os.path.join(results_dir, f"{base_filename}_{ff_type}.data")  # Use .data extension
                # Load forcefield parameters (Pair Coeffs) from JSON for LAMMPS output.
                # Select the mineral block matching the chosen k_angle and forcefield variant.
                try:
                    if ff_type == 'clayff':
                        mineral_block = 'CLAYFF_2004'
                    else:
                        k_map = {'0': 'GMINFF_k0', '250': 'GMINFF_k250', '500': 'GMINFF_k500', '1500': 'GMINFF_k1500'}
                        mineral_block = k_map.get(str(angle_terms), 'GMINFF_k500')
                    ap = get_ap()
                    ff_params = ap.load_forcefield(
                        'GMINFF/gminff_all.json',
                        blocks=[mineral_block, 'OPC3_HFE_LM', 'OPC3']
                    )
                except Exception as e:
                    print(f"Warning: Could not load forcefield Pair Coeffs: {e}")
                    ff_params = None
                from atomipy import write_top
                write_top.lmp(
                    atoms,
                    Box=Box_dim,
                    file_path=topology_lmp,
                    forcefield=ff_params,
                    KANGLE=topology_kangle,
                    max_angle=topology_max_angle,
                )
                generated_files.append(topology_lmp)
                print(f"LAMMPS data file for {ff_type} written successfully to {topology_lmp}")
                progress_step += progress_increment
        else:
            # If topology generation is disabled, skip to the next step and update progress
            tasks_status[task_id] = {'status': 'Processing', 'step': 'Skipping topology files (disabled)', 'progress': 85}
            print("Topology generation skipped as requested.")
        # PDB file is always generated for reference
        # If topology was disabled, progress is already at 85, so don't update it again
        if generate_topology or progress_step < 85:
            tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing PDB', 'progress': 85}
        
        pdb_filepath = os.path.join(results_dir, f"{base_filename}_{ff_type}.pdb")
        # Convert Box_dim to Cell parameters for PDB and XYZ formats
        ap = get_ap()
        Cell = ap.Box_dim2Cell(Box_dim)
        ap.write_pdb(atoms, Cell, pdb_filepath)
        print(f"PDB file written successfully to {pdb_filepath}")

        # GRO file is also always generated for reference
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing GRO', 'progress': 90} # Update progress
        gro_filepath = os.path.join(results_dir, f"{base_filename}_{ff_type}.gro")
        ap = get_ap()
        ap.write_gro(atoms, Box=Box_dim, file_path=gro_filepath)
        print(f"GRO file written successfully to {gro_filepath}")
        
        # XYZ file is also generated for reference with Cell info on line 2
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing XYZ', 'progress': 95} # Final writing step
        xyz_filepath = os.path.join(results_dir, f"{base_filename}_{ff_type}.xyz")
        ap = get_ap()
        ap.write_xyz(atoms, Box=Cell, file_path=xyz_filepath)
        print(f"XYZ file written successfully to {xyz_filepath}")

        elapsed_time = time.time() - start_time
        logger.info(f"Task {task_id} completed successfully in {elapsed_time:.2f} seconds. Results ID: {results_id}")
        tasks_status[task_id] = {
            'status': 'Complete', 
            'results_id': results_id, 
            'progress': 100,
            'elapsed_time': f"{elapsed_time:.2f} seconds"
        }
    except Exception as e:
        error_msg = f"Error processing task {task_id}: {str(e)}"
        logger.error(error_msg)
        logger.error(traceback.format_exc())  # Log detailed traceback for debugging
        
        # Try to extract more details about the error
        error_type = type(e).__name__
        error_details = str(e)
        
        # Create a more informative error message
        detailed_error = f"{error_type}: {error_details}"
        logger.error(f"Detailed error: {detailed_error}")
        
        tasks_status[task_id] = {
            'status': 'Error', 
            'message': f'Error generating output files: {detailed_error}', 
            'progress': 100,
            'error_type': error_type,
            'timestamp': time.time()
        }
    finally:
        # Clean up the original uploaded file from the results directory
        # as it's not needed after processing.
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except OSError as e:
                print(f"Error removing uploaded file {filepath}: {e}")

# --- Routes ---
@app.route('/')
def index():
    # Clean up old sessions? Or maybe rely on temp dir cleanup
    return render_template('index.html', gemmi_available=app.config['GEMMI_AVAILABLE'])

@app.route('/upload_file', methods=['POST'])
def start_processing_task():  # Renamed route function
    if 'file' not in request.files:
        flash('No file part')
        return jsonify({'error': 'No file part'}), 400  # Return JSON error
    file = request.files['file']
    if file.filename == '':
        flash('No selected file')
        return jsonify({'error': 'No selected file'}), 400  # Return JSON error

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        base_filename = filename.rsplit('.', 1)[0]
        file_extension = filename.rsplit('.', 1)[1].lower()

        if file_extension in {'cif', 'mmcif', 'mcif'} and not app.config['GEMMI_AVAILABLE']:
            return jsonify({
                'error': (
                    'CIF support requires GEMMI, but GEMMI is not installed on this server. '
                    'Install with: pip install gemmi'
                )
            }), 400

        # Create a unique results directory for this task
        results_id = str(uuid.uuid4())
        results_dir = os.path.join(app.config['RESULTS_FOLDER'], results_id)
        os.makedirs(results_dir, exist_ok=True)

        # Save the uploaded file temporarily (could save directly to results?)
        # filepath = os.path.join(session_dir, filename)
        filepath = os.path.join(results_dir, filename)  # Save in results dir
        file.save(filepath)

        ff_type = request.form.get('forcefield', 'minff')
        # Symmetry expansion is always enabled for CIF import in the web app.
        expand_symmetry = True
        fuse_overlaps = parse_bool(request.form.get('fuse_overlaps'), default=False)
        fuse_rmax = parse_float(request.form.get('fuse_rmax'), DEFAULT_CIF_FUSE_RMAX)
        fuse_criteria = request.form.get('fuse_criteria', DEFAULT_CIF_FUSE_CRITERIA).strip().lower()
        add_hydrogen_bvs = parse_bool(request.form.get('add_hydrogen_bvs'), default=False)
        bvs_delta_threshold = parse_float(
            request.form.get('bvs_delta_threshold'),
            DEFAULT_BVS_DELTA_THRESHOLD,
        )
        bvs_max_additions = parse_int(
            request.form.get('bvs_max_additions'),
            DEFAULT_BVS_MAX_ADDITIONS,
        )
        bvs_max_additions = max(0, bvs_max_additions)
        replicate_structure = parse_bool(request.form.get('replicate_structure'), default=False)
        replicate_nx = max(1, parse_int(request.form.get('replicate_nx'), 1))
        replicate_ny = max(1, parse_int(request.form.get('replicate_ny'), 1))
        replicate_nz = max(1, parse_int(request.form.get('replicate_nz'), 1))
        reset_molid = parse_bool(request.form.get('reset_molid'), default=True)
        default_angle_terms = DEFAULT_ANGLE_TERM_MINFF if ff_type == 'minff' else DEFAULT_ANGLE_TERM_CLAYFF
        angle_terms = request.form.get('angle_terms', default_angle_terms).strip().lower()
        if angle_terms not in ANGLE_TERM_OPTIONS:
            angle_terms = default_angle_terms

        if fuse_criteria not in {'average', 'occupancy', 'order'}:
            fuse_criteria = DEFAULT_CIF_FUSE_CRITERIA
        if fuse_rmax <= 0:
            fuse_rmax = DEFAULT_CIF_FUSE_RMAX
        
        # Check if topology generation is enabled
        generate_topology = request.form.get('generate_topology', 'true').lower() == 'true'
        
        # Only validate output formats if topology generation is enabled
        output_formats = []
        if generate_topology:
            output_formats = request.form.getlist('output_formats')
            if not output_formats:
                # Clean up created results dir if no format selected when topology is enabled
                if os.path.exists(results_dir):
                    try:
                        shutil.rmtree(results_dir)
                    except OSError as e:
                        print(f"Error removing directory {results_dir}: {e}")
                return jsonify({'error': 'Please select at least one output topology format or disable topology generation.'}), 400

        try:
            # Generate a unique task ID
            task_id = str(uuid.uuid4())
            tasks_status[task_id] = {'status': 'Pending', 'progress': 0}

            # Submit the task to the executor with a timeout
            try:
                # Log the task submission
                file_size = os.path.getsize(filepath)
                logger.info(f"Submitting task {task_id} for file {filename} ({file_size} bytes)")
                
                # Submit the task to the executor
                executor.submit(
                    process_file_task,
                    task_id,
                    filepath,
                    filename,
                    ff_type,
                    output_formats,
                    results_id,
                    results_dir,
                    generate_topology,
                    expand_symmetry,
                    fuse_overlaps,
                    fuse_rmax,
                    fuse_criteria,
                    add_hydrogen_bvs,
                    bvs_delta_threshold,
                    bvs_max_additions,
                    replicate_structure,
                    replicate_nx,
                    replicate_ny,
                    replicate_nz,
                    angle_terms,
                    reset_molid,
                )
            except Exception as e:
                logger.error(f"Error submitting task: {str(e)}")
                raise

            # Return the task ID to the client
            return jsonify({'task_id': task_id})
        except RequestEntityTooLarge:
            flash('File too large. Maximum size allowed is 16 MB.')
            return jsonify({'error': 'File too large'}), 413
        except Exception as e:  # Catch potential errors before task submission (e.g., saving file)
            print(f"Error preparing task: {e}")
            # Clean up results dir if created
            if os.path.exists(results_dir):
                try:
                    shutil.rmtree(results_dir)
                except OSError as e:
                    print(f"Error removing directory {results_dir}: {e}")
            return jsonify({'error': f'An error occurred preparing the task: {e}'}), 500

    else:
        flash('Invalid file type. Allowed types are .gro, .pdb, .xyz, .cif')
        return jsonify({'error': 'Invalid file type'}), 400

@app.route('/results/<results_id>')
def results(results_id):
    results_dir = os.path.join(app.config['RESULTS_FOLDER'], results_id)
    if not os.path.exists(results_dir):
        flash(f'Results not found for ID {results_id} or may have expired.')
        return redirect(url_for('index'))  # Redirect to index if results don't exist

    # Get list of all result files
    try:
        result_files = os.listdir(results_dir)
        result_files = [f for f in result_files if os.path.isfile(os.path.join(results_dir, f))]
    except Exception as e:
        flash(f'Could not find results directory: {e}')
        return redirect(url_for('index'))

    # Organize files by type and forcefield
    minff_files = [f for f in result_files if '_minff.' in f]
    clayff_files = [f for f in result_files if '_clayff.' in f]
    
    # Handle log files explicitly
    minff_logs = [f for f in result_files if 'minff_structure_stats.log' in f]
    clayff_logs = [f for f in result_files if 'clayff_structure_stats.log' in f]
    
    # Any file not caught by the above categories
    other_files = [f for f in result_files if 
                  ('_minff.' not in f and '_clayff.' not in f) and
                  ('minff_structure_stats.log' not in f and 'clayff_structure_stats.log' not in f)]

    files = {
        'minff_structures': [f for f in minff_files if f.endswith('.gro') or f.endswith('.pdb') or f.endswith('.xyz')],
        'minff_topologies': [f for f in minff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
        'minff_logs': minff_logs,
        'clayff_structures': [f for f in clayff_files if f.endswith('.gro') or f.endswith('.pdb') or f.endswith('.xyz')],
        'clayff_topologies': [f for f in clayff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
        'clayff_logs': clayff_logs,
        'other_files': other_files
    }

    return render_template('results.html', results_id=results_id, files=files)

@app.route('/status/<task_id>')
def get_status(task_id):
    status = tasks_status.get(task_id, {'status': 'Unknown', 'progress': 0, 'message': 'Task ID not found.'})
    return jsonify(status)

@app.route('/task_result/<task_id>')
def get_task_result(task_id):
    status = tasks_status.get(task_id)
    if status and status.get('status') == 'Complete':
        results_id = status.get('results_id')
        if results_id:
            # Optionally clear task from memory after retrieval?
            # tasks_status.pop(task_id, None)
            return redirect(url_for('results', results_id=results_id))
        else:
            flash('Task completed but results ID is missing.')
            return redirect(url_for('index'))
    elif status and status.get('status') == 'Error':
        flash(f"Processing failed: {status.get('message', 'Unknown error')}")
        # Optionally clear task from memory after retrieval?
        # tasks_status.pop(task_id, None)
        return redirect(url_for('index'))
    elif status:
        # Task is still processing or in an unknown state - should not redirect here from button click
        # This might happen if user manually navigates here
        flash(f'Task {task_id} is still processing or in an unknown state.')
        return redirect(url_for('index'))  # Or maybe a dedicated 'processing' page?
    else:
        flash(f'Unknown Task ID: {task_id}')
        return redirect(url_for('index'))

@app.route('/download/<results_id>/<filename>')
def download_file(results_id, filename):
    directory = os.path.join(app.config['RESULTS_FOLDER'], results_id)
    return send_from_directory(directory, filename, as_attachment=True)

@app.route('/download_zip/<results_id>')
def download_zip(results_id):
    results_dir = os.path.join(app.config['RESULTS_FOLDER'], results_id)
    if not os.path.exists(results_dir):
        flash('Results not found.')
        return redirect(url_for('index'))

    # Create a list of all files in the results directory
    file_list = [f for f in os.listdir(results_dir) if os.path.isfile(os.path.join(results_dir, f))]

    # Create a ZIP file in memory
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w') as zf:
        for file in file_list:
            file_path = os.path.join(results_dir, file)
            zf.write(file_path, arcname=file)

    memory_file.seek(0)

    # Create a unique filename for the download
    zip_filename = f"atomipy_results_{results_id[:8]}.zip"

    # Return the ZIP file as a response
    return Response(
        memory_file,
        mimetype="application/zip",
        headers={
            "Content-Disposition": f"attachment;filename={zip_filename}"
        }
    )

@app.route('/about')
def about():
    return render_template('about.html', gemmi_available=app.config['GEMMI_AVAILABLE'])

import threading
import contextlib
import shutil

def prune_cache_loop():
    """Background thread to delete result files older than 1 hour."""
    while True:
        try:
            now = time.time()
            cutoff = now - 3600 # 1 hour
            # We access app.config through a manual folder check since it's a global app
            folders = [os.path.join(BASE_DIR, 'uploads'), os.path.join(BASE_DIR, 'results')]
            for folder in folders:
                if os.path.exists(folder):
                    for f in os.listdir(folder):
                        p = os.path.join(folder, f)
                        if os.path.getmtime(p) < cutoff:
                            with contextlib.suppress(Exception):
                                if os.path.isfile(p):
                                    os.remove(p)
                                elif os.path.isdir(p):
                                    shutil.rmtree(p)
        except Exception as e:
            print(f"Error in pruning thread: {e}")
        time.sleep(1800) # Run every 30 mins

# Start the pruning thread
threading.Thread(target=prune_cache_loop, daemon=True).start()

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=5001)
