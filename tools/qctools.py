# Tools for molecular modelling with Jupyter Notebooks
# 2021-2026
# antti.karttunen@aalto.fi

############## Functions to facilitate printing ##############

def print_info(info):
    # info: string
    banner = "-----------------------------------------------------------"
    print(f"{banner}\n{info}\n{banner}")
          
def print_error(error):
    # error: string
    err = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"      
    print(f"{err}\n{error}\n{err}")

############## Functions for loading molecules ##############

def load_xyz(xyzfile, silent = False):
    """Loads XYZ file
    xyzfile: file name or path as string
    silent: no printing
    
    Returns: ase.Atoms or None if fails
    """
    from ase.io import read    
    try:
        atoms = read(xyzfile, format = 'xyz')
        if not silent:
            print_info(f"Molecule was loaded from file {xyzfile}\n"
                       f"Atoms: {len(atoms)}\n"
                       f"Formula: {str(atoms.symbols)}")
    except OSError:
        if not silent:
            print_error(f"Failed to load file {xyzfile}")
        return None
    else:
        atoms.info['id'] = xyzfile
        return atoms

def load_xyz_as_traj(xyztraj, silent = False):
    """ Converts multi-XYZ file to ASE trajectory (to visualize with nglview)
    xyztraj: file name or path as a string (multi-XYZ file)    
    silent: no printing

    Returns: ase.io.Trajectory or None if fails
    """
    from ase.io import read, write, Trajectory

    # Read the XYZ file (index=':' reads all frames/images)
    try:
        frames = read(xyztraj, index=':')
        if not silent:
            print_info(f"XYZ trajectory was loaded from file {xyztraj}\n"
                       f"Frames: {len(frames)}")
    except OSError:
        if not silent:
            print_error(f"Failed to load file {xyzfile}")
        return None        
    
    # Write the frames to a temporary .traj file and load as Trajectory
    try:
        trajfile = xyztraj + '.traj'
        with Trajectory(trajfile, mode='w') as traj:
            for frame in frames:
                traj.write(frame)            
    except OSError:
        if not silent:
            print_error(f"Failed to save file {trajfile}")
        return None
        
    try:
        traj = Trajectory(trajfile)
    except OSError:
        if not silent:
            print_error(f"Failed to load trajectory file {trajfile}")
    else:
        return traj
            
############## Functions for visualizing molecules with nglview

def show_molecule(molecule, size = (500, 400), style = 'ball+stick', unitcell = None, labels = None, bg = 'black'):
    """ Shows a molecule using NGLWidget
    molecule: ase.Atoms or ase.io.Trajectory
    size: tuple of two integers (x_pixels, y_pixels)
    style: 'ball+stick', 'spacefill', 'licorice', 'line' (http://nglviewer.org/ngl/api/manual/molecular-representations.html)
    unitcell: None or 'white', 'orange', 'red', ...
    labels: None or 'atomname', 'atomindex, 'element'
    bg: 'black', 'white', ...
    Returns: NGLWidget
    """
    import nglview
    from ase import Atoms
    from ase.io.trajectory import TrajectoryReader

    if isinstance(molecule, Atoms):
        nv = nglview.show_ase(molecule)
    elif isinstance(molecule, TrajectoryReader):
        nv = nglview.show_asetraj(molecule)
    else:
        print_error("Invalid molecule!")
        return None
    
    nv._set_size(f"{size[0]}px", f"{size[1]}px")
    # Calls: nv._remote_call('setSize', target='Widget', args=[w, h])
    
    nv.clear_representations()    
    nv.add_representation(style)
    
    if style == 'spacefill':
        # radiusType = 'covalent' looks better for bulk materials. Can be changed to 'vdw' via GUI
        nv.update_representation(component = 0, repr_index = 0, radiusType='covalent')

    if unitcell is not None:
        nv.add_representation('unitcell')        
        # Changing unit cell color does not work, not even via NGL GUI?
        # nv.add_representation('unitcell', colorValue = '#ffffff')
    if labels is not None:
        nv.add_representation('label', labelType = labels)
        
    nv.parameters = dict(backgroundColor = bg, clipDist = -100)
    nv.camera = 'orthographic'
    nv.display()
    return nv  

