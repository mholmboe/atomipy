import os
import sys
import uuid
import zipfile
import io
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash, session, Response, jsonify
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge
from flask_executor import Executor
import shutil

# Add the parent directory to sys.path to import atomipy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import atomipy as a module for access to all functionality
import atomipy as ap

# Direct imports for functions explicitly exposed in __init__.py
from atomipy import import_gro, import_pdb
from atomipy import write_conf, write_top
from atomipy import Box_dim2Cell, Cell2Box_dim

# Updated imports for forcefield and charge modules
from atomipy.forcefield import minff, clayff

app = Flask(__name__)
app.config['SECRET_KEY'] = os.urandom(24)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB upload limit

# Initialize Flask-Executor
executor = Executor(app)
tasks_status = {}  # Dictionary to store task status

# Define upload and results folders
UPLOAD_FOLDER_NAME = 'uploads'
RESULTS_FOLDER_NAME = 'results'

# Use /tmp for writable storage on App Engine, otherwise use local folders
if os.environ.get('GAE_ENV', '').startswith('standard'):
    app.config['UPLOAD_FOLDER'] = '/tmp/uploads'
    app.config['RESULTS_FOLDER'] = '/tmp/results'
else:
    app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), UPLOAD_FOLDER_NAME)
    app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), RESULTS_FOLDER_NAME)

app.config['ALLOWED_EXTENSIONS'] = {'gro', 'pdb'}

# Create upload and results folders if they don't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

# --- Background Task ---
def process_file_task(task_id, filepath, filename, ff_type, output_formats, results_id, results_dir):
    base_filename = filename.rsplit('.', 1)[0]
    try:
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Reading structure', 'progress': 10}
        # Use appropriate import function based on file extension
        file_extension = filename.rsplit('.', 1)[1].lower()
        if file_extension == 'gro':
            atoms, box_dim = import_gro(filepath)
        elif file_extension == 'pdb':
            atoms, cell = import_pdb(filepath)
            # Convert cell to box_dim for consistency
            box_dim = Cell2Box_dim(cell)
        else:
            raise ValueError(f"Unsupported file extension: {file_extension}")
        tasks_status[task_id] = {'status': 'Processing', 'step': f'Assigning {ff_type} atom types', 'progress': 30}
        if ff_type == 'minff':
            atoms = minff(atoms, box_dim)
        elif ff_type == 'clayff':
            atoms = clayff(atoms, box_dim)
        tasks_status[task_id] = {'status': 'Processing', 'step': f'Calculating charges ({ff_type})', 'progress': 50}
        # Comprehensive debug of the structure
        print(f"Type of atoms after processing: {type(atoms)}")
        
        # Special handling for tuple returned from forcefield functions
        if isinstance(atoms, tuple) and len(atoms) == 2:
            # If atoms is a tuple of (atoms_list, box_dim), extract the atoms list
            atoms_from_tuple, box_dim_from_tuple = atoms
            atoms = atoms_from_tuple
            print(f"Extracted atoms list from tuple, new box_dim: {box_dim_from_tuple}")
            # Update box_dim if needed
            if box_dim_from_tuple is not None and len(box_dim_from_tuple) in [3, 9]:
                box_dim = box_dim_from_tuple
                print(f"Updated box_dim to: {box_dim}")
        
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
        # Generate selected output files
        generated_files = []
        progress_step = 70
        progress_increment = (95 - progress_step) / len(output_formats) if output_formats else 0
        if 'itp' in output_formats:
            tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing ITP', 'progress': int(progress_step)}
            topology_itp = os.path.join(results_dir, f"{base_filename}_{ff_type}.itp")
            write_top.itp(atoms, box_dim, file_path=topology_itp)
            generated_files.append(topology_itp)
            print(f"ITP file for {ff_type} written successfully to {topology_itp}")
            progress_step += progress_increment
        if 'psf' in output_formats:
            tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing PSF', 'progress': int(progress_step)}
            topology_psf = os.path.join(results_dir, f"{base_filename}_{ff_type}.psf")
            write_top.psf(atoms, box_dim, file_path=topology_psf)
            generated_files.append(topology_psf)
            print(f"PSF file for {ff_type} written successfully to {topology_psf}")
            progress_step += progress_increment
        if 'lmp' in output_formats:  # Check for 'lmp' from form value
            tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing LAMMPS Data', 'progress': int(progress_step)}
            topology_lmp = os.path.join(results_dir, f"{base_filename}_{ff_type}.data")  # Use .data extension
            write_top.lmp(atoms, box_dim, file_path=topology_lmp)
            generated_files.append(topology_lmp)
            print(f"LAMMPS data file for {ff_type} written successfully to {topology_lmp}")
            progress_step += progress_increment
        # PDB file is always generated for reference
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing PDB', 'progress': 90} # Adjust progress slightly
        pdb_filepath = os.path.join(results_dir, f"{base_filename}_{ff_type}.pdb")
        write_conf.pdb(atoms, Box_dim2Cell(box_dim), pdb_filepath)
        print(f"PDB file written successfully to {pdb_filepath}")

        # GRO file is also always generated for reference
        tasks_status[task_id] = {'status': 'Processing', 'step': 'Writing GRO', 'progress': 95} # Final writing step
        gro_filepath = os.path.join(results_dir, f"{base_filename}_{ff_type}.gro")
        write_conf.gro(atoms, box_dim, gro_filepath)
        print(f"GRO file written successfully to {gro_filepath}")

        tasks_status[task_id] = {'status': 'Complete', 'results_id': results_id, 'progress': 100}
        print(f"Task {task_id} completed successfully. Results ID: {results_id}")
    except Exception as e:
        print(f"Error processing task {task_id}: {e}")
        import traceback
        traceback.print_exc()  # Log traceback for debugging
        tasks_status[task_id] = {'status': 'Error', 'message': f'Error generating output files: {e}', 'progress': 100}
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
    return render_template('index.html')

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

        # Create a unique results directory for this task
        results_id = str(uuid.uuid4())
        results_dir = os.path.join(app.config['RESULTS_FOLDER'], results_id)
        os.makedirs(results_dir, exist_ok=True)

        # Save the uploaded file temporarily (could save directly to results?)
        # filepath = os.path.join(session_dir, filename)
        filepath = os.path.join(results_dir, filename)  # Save in results dir
        file.save(filepath)

        ff_type = request.form.get('forcefield', 'minff')
        output_formats = request.form.getlist('output_formats')

        if not output_formats:
            # Clean up created results dir if no format selected
            if os.path.exists(results_dir):
                try:
                    shutil.rmtree(results_dir)
                except OSError as e:
                    print(f"Error removing directory {results_dir}: {e}")
            return jsonify({'error': 'Please select at least one output topology format.'}), 400

        try:
            # Generate a unique task ID
            task_id = str(uuid.uuid4())
            tasks_status[task_id] = {'status': 'Pending', 'progress': 0}

            # Submit the task to the executor
            executor.submit(process_file_task, task_id, filepath, filename, ff_type, output_formats, results_id, results_dir)

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
        flash('Invalid file type. Allowed types are .gro and .pdb')
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
    other_files = [f for f in result_files if '_minff.' not in f and '_clayff.' not in f]

    files = {
        'minff_structures': [f for f in minff_files if f.endswith('.gro') or f.endswith('.pdb')],
        'minff_topologies': [f for f in minff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
        'clayff_structures': [f for f in clayff_files if f.endswith('.gro') or f.endswith('.pdb')],
        'clayff_topologies': [f for f in clayff_files if f.endswith('.itp') or f.endswith('.psf') or f.endswith('.data')],
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
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)
