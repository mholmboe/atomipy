import os
import sys
import uuid
import zipfile
import io
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash, session, Response
from werkzeug.utils import secure_filename

# Add the parent directory to sys.path to import atomipy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import atomipy as a module for access to all functionality
import atomipy as ap

# Direct imports for functions explicitly exposed in __init__.py
from atomipy import import_gro, import_pdb
from atomipy import write_conf, write_top
from atomipy import element, bond_angle
from atomipy import Box_dim2Cell, Cell2Box_dim

# Updated imports for forcefield and charge modules
from atomipy.forcefield import minff, clayff
from atomipy.charge import charge_minff, charge_clayff, assign_formal_charges, balance_charges

app = Flask(__name__)
app.config['SECRET_KEY'] = os.urandom(24)

# Use /tmp for writable storage on App Engine, otherwise use local folders
if os.environ.get('GAE_ENV', '').startswith('standard'):
    app.config['UPLOAD_FOLDER'] = '/tmp/uploads'
    app.config['RESULTS_FOLDER'] = '/tmp/results'
else:
    app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uploads')
    app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')

app.config['ALLOWED_EXTENSIONS'] = {'gro', 'pdb'}
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB max upload size

# Create upload and results folders if they don't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        flash('No file part')
        return redirect(request.url)
    
    file = request.files['file']
    if file.filename == '':
        flash('No selected file')
        return redirect(request.url)
    
    if file and allowed_file(file.filename):
        # Generate a unique identifier for this session
        session_id = str(uuid.uuid4())
        session['session_id'] = session_id
        
        # Create a session directory for this upload
        session_dir = os.path.join(app.config['UPLOAD_FOLDER'], session_id)
        results_dir = os.path.join(app.config['RESULTS_FOLDER'], session_id)
        os.makedirs(session_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)
        
        # Save the uploaded file
        filename = secure_filename(file.filename)
        file_extension = filename.rsplit('.', 1)[1].lower()
        filepath = os.path.join(session_dir, filename)
        file.save(filepath)
        
        # Get forcefield option from the form
        ff_type = request.form.get('forcefield', 'minff')
        
        try:
            # Process the file based on its extension
            try:
                if file_extension == 'gro':
                    atoms, box_dim = import_gro(filepath)
                elif file_extension == 'pdb':
                    atoms, cell = import_pdb(filepath)
                    # Convert cell to box_dim for consistency
                    box_dim = Cell2Box_dim(cell)
                else:
                    flash('Unsupported file format')
                    return redirect(url_for('index'))
                
                # Print debugging info
                print(f"Loaded {len(atoms)} atoms from {filename}")
                print(f"Sample atom structure: {atoms[0] if atoms else 'No atoms'}")
                print(f"Box dimensions: {box_dim}")
                print(f"Selected forcefield: {ff_type}")
                print(f"Charges will be assigned automatically based on forcefield: {ff_type}")
                
            except Exception as e:
                flash(f'Error importing file: {str(e)}')
                import traceback
                traceback.print_exc()
                return redirect(url_for('index'))
            
            # Perform preprocessing - convert tuples to dictionaries if needed
            try:
                print("Converting any tuple atoms to dictionaries...")
                for i in range(len(atoms)):
                    if not isinstance(atoms[i], dict):
                        print(f"Atom {i} is not a dictionary, but a {type(atoms[i])}")
                        # If it's a tuple, convert to dictionary
                        if isinstance(atoms[i], tuple):
                            print(f"Converting tuple atom {i} to dictionary")
                            atom_dict = {}
                            # Typical atom properties
                            properties = ['aid', 'type', 'resid', 'resname', 'aname', 'q', 'mass', 'x', 'y', 'z']
                            for j, prop in enumerate(properties):
                                if j < len(atoms[i]):
                                    atom_dict[prop] = atoms[i][j]
                            atoms[i] = atom_dict
                        else:
                            raise TypeError(f"Cannot process atom {i}: expected dict, got {type(atoms[i])}")
                print("Preprocessing completed successfully")
            except Exception as e:
                flash(f'Error during preprocessing: {str(e)}')
                import traceback
                traceback.print_exc()
                return redirect(url_for('index'))
            
            # Apply forcefield based on user selection
            try:
                if ff_type == 'minff':
                    print("Starting MINFF atom type assignment...")
                    # Note: minff() can return either just atoms or a tuple of (atoms, box_dim)
                    minff_result = minff(atoms, box_dim)
                    
                    # Check if minff returned a tuple (atoms, box_dim) or just atoms
                    if isinstance(minff_result, tuple) and len(minff_result) == 2:
                        print("minff() returned a tuple of (atoms, box_dim)")
                        atoms, box_dim_updated = minff_result
                        print(f"Updated box_dim: {box_dim_updated}")
                    else:
                        print("minff() returned just the atoms list")
                        atoms = minff_result
                        
                    print("MINFF atom type assignment completed successfully")
                elif ff_type == 'clayff':
                    print("Starting CLAYFF atom type assignment...")
                    # Note: clayff() can return either just atoms or a tuple of (atoms, box_dim)
                    clayff_result = clayff(atoms, box_dim)
                    
                    # Check if clayff returned a tuple (atoms, box_dim) or just atoms
                    if isinstance(clayff_result, tuple) and len(clayff_result) == 2:
                        print("clayff() returned a tuple of (atoms, box_dim)")
                        atoms, box_dim_updated = clayff_result
                        print(f"Updated box_dim: {box_dim_updated}")
                    else:
                        print("clayff() returned just the atoms list")
                        atoms = clayff_result
                        
                    print("CLAYFF atom type assignment completed successfully")
                else:
                    flash('Invalid forcefield type selected')
                    return redirect(url_for('index'))
                
                # Comprehensive debug of the structure
                print(f"Type of atoms after processing: {type(atoms)}")
                
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
            except Exception as e:
                flash(f'Error during forcefield processing: {str(e)}')
                import traceback
                traceback.print_exc()
                return redirect(url_for('index'))
                        # 4. Generate output files with new atom types
            try:
                print("Starting to generate output files...")
                base_filename = os.path.splitext(filename)[0]
                
                # Save processed GRO file for primary forcefield
                print(f"Writing GRO file for {ff_type}...")
                processed_gro = os.path.join(results_dir, f"{base_filename}_{ff_type}.gro")
                write_conf.gro(atoms, box_dim, processed_gro)
                print(f"GRO file written successfully to {processed_gro}")
                
                # Save processed PDB file for primary forcefield
                print(f"Writing PDB file for {ff_type}...")
                processed_pdb = os.path.join(results_dir, f"{base_filename}_{ff_type}.pdb")
                # Convert box_dim back to cell for PDB
                cell = Box_dim2Cell(box_dim)
                write_conf.pdb(atoms, cell, processed_pdb)
                print(f"PDB file written successfully to {processed_pdb}")
                
                # Create a deep copy of the original atoms for the secondary forcefield
                import copy
                original_atoms = None
                
                # Store the current processed atoms based on the selected forcefield
                if file_extension == 'gro':
                    # Re-import the original file to get a clean slate
                    original_atoms, original_box_dim = import_gro(filepath)
                elif file_extension == 'pdb':
                    # Re-import the original file to get a clean slate
                    original_atoms, original_cell = import_pdb(filepath)
                    original_box_dim = Cell2Box_dim(original_cell)
                
                # Apply the other forcefield to generate alternative topology files
                alt_ff_type = 'clayff' if ff_type == 'minff' else 'minff'
                
                # Process with the alternative forcefield
                alt_atoms = copy.deepcopy(original_atoms)
                print(f"Processing with alternative forcefield: {alt_ff_type}")
                
                try:
                    if alt_ff_type == 'minff':
                        alt_result = minff(alt_atoms, original_box_dim)
                        if isinstance(alt_result, tuple) and len(alt_result) == 2:
                            alt_atoms, alt_box_dim = alt_result
                        else:
                            alt_atoms = alt_result
                            alt_box_dim = original_box_dim
                    else:  # alt_ff_type == 'clayff'
                        alt_result = clayff(alt_atoms, original_box_dim)
                        if isinstance(alt_result, tuple) and len(alt_result) == 2:
                            alt_atoms, alt_box_dim = alt_result
                        else:
                            alt_atoms = alt_result
                            alt_box_dim = original_box_dim
                    
                    # Save structure files for alternative forcefield
                    print(f"Writing GRO file for {alt_ff_type}...")
                    alt_processed_gro = os.path.join(results_dir, f"{base_filename}_{alt_ff_type}.gro")
                    write_conf.gro(alt_atoms, alt_box_dim, alt_processed_gro)
                    
                    print(f"Writing PDB file for {alt_ff_type}...")
                    alt_cell = Box_dim2Cell(alt_box_dim)
                    alt_processed_pdb = os.path.join(results_dir, f"{base_filename}_{alt_ff_type}.pdb")
                    write_conf.pdb(alt_atoms, alt_cell, alt_processed_pdb)
                except Exception as e:
                    print(f"Warning: Could not process with alternative forcefield {alt_ff_type}: {str(e)}")
                    alt_atoms = None
                
                # Generate topology files for primary forcefield
                # GROMACS (.itp)
                print(f"Writing GROMACS topology file for {ff_type}...")
                topology_itp = os.path.join(results_dir, f"{base_filename}_{ff_type}.itp")
                write_top.itp(atoms, box_dim, topology_itp)
                print(f"ITP file for {ff_type} written successfully to {topology_itp}")
                
                # NAMD (.psf)
                print(f"Writing NAMD topology file for {ff_type}...")
                topology_psf = os.path.join(results_dir, f"{base_filename}_{ff_type}.psf")
                write_top.psf(atoms, box_dim, topology_psf)
                print(f"PSF file for {ff_type} written successfully to {topology_psf}")
                
                # LAMMPS (.data)
                print(f"Writing LAMMPS data file for {ff_type}...")
                topology_data = os.path.join(results_dir, f"{base_filename}_{ff_type}.data")
                write_top.lmp(atoms, box_dim, topology_data)
                print(f"LAMMPS data file for {ff_type} written successfully to {topology_data}")
                
                # Generate topology files for alternative forcefield if processing was successful
                if alt_atoms:
                    # GROMACS (.itp)
                    print(f"Writing GROMACS topology file for {alt_ff_type}...")
                    alt_topology_itp = os.path.join(results_dir, f"{base_filename}_{alt_ff_type}.itp")
                    write_top.itp(alt_atoms, alt_box_dim, alt_topology_itp)
                    print(f"ITP file for {alt_ff_type} written successfully to {alt_topology_itp}")
                    
                    # NAMD (.psf)
                    print(f"Writing NAMD topology file for {alt_ff_type}...")
                    alt_topology_psf = os.path.join(results_dir, f"{base_filename}_{alt_ff_type}.psf")
                    write_top.psf(alt_atoms, alt_box_dim, alt_topology_psf)
                    print(f"PSF file for {alt_ff_type} written successfully to {alt_topology_psf}")
                    
                    # LAMMPS (.data)
                    print(f"Writing LAMMPS data file for {alt_ff_type}...")
                    alt_topology_data = os.path.join(results_dir, f"{base_filename}_{alt_ff_type}.data")
                    write_top.lmp(alt_atoms, alt_box_dim, alt_topology_data)
                    print(f"LAMMPS data file for {alt_ff_type} written successfully to {alt_topology_data}")
                
                print("All output files generated successfully")
            except Exception as e:
                flash(f'Error generating output files: {str(e)}')
                import traceback
                traceback.print_exc()
                return redirect(url_for('index'))
            
            # Return the results page with links to download the files
            return redirect(url_for('results', session_id=session_id))
            
        except Exception as e:
            flash(f'Error processing file: {str(e)}')
            return redirect(url_for('index'))
    
    flash('Invalid file format. Please upload a .gro or .pdb file.')
    return redirect(url_for('index'))

@app.route('/results/<session_id>')
def results(session_id):
    results_dir = os.path.join(app.config['RESULTS_FOLDER'], session_id)
    if not os.path.exists(results_dir):
        flash('Results not found.')
        return redirect(url_for('index'))
    
    # Get list of all result files
    result_files = os.listdir(results_dir)
    result_files = [f for f in result_files if os.path.isfile(os.path.join(results_dir, f))]
    
    # Organize files by type and forcefield
    minff_files = [f for f in result_files if '_minff.' in f]
    clayff_files = [f for f in result_files if '_clayff.' in f]
    other_files = [f for f in result_files if '_minff.' not in f and '_clayff.' not in f]
    
    files = {
        'minff_structures': [f for f in minff_files if f.endswith('.gro') or f.endswith('.pdb')],
        'minff_topologies': [f for f in minff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
        'clayff_structures': [f for f in clayff_files if f.endswith('.gro') or f.endswith('.pdb')],
        'clayff_topologies': [f for f in clayff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
        'other_files': other_files
    }
    
    return render_template('results.html', session_id=session_id, files=files)

@app.route('/download/<session_id>/<filename>')
def download_file(session_id, filename):
    results_dir = os.path.join(app.config['RESULTS_FOLDER'], session_id)
    return send_from_directory(results_dir, filename, as_attachment=True)

@app.route('/download_zip/<session_id>')
def download_zip(session_id):
    results_dir = os.path.join(app.config['RESULTS_FOLDER'], session_id)
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
    zip_filename = f"atomipy_results_{session_id[:8]}.zip"
    
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
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)
