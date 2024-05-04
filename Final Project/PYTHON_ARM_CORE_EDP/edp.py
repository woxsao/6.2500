# ARM Core energy and delay simulator for 6.2500 class project.
# 
# Class `ARMCore` is main simulation class. It loads ARM core layout data
# from hdf5 file and calculates energy and delay of the core using
# configurable MOSFET parameters found from Sentaurus TCAD simulation.
#
# Class `EDPGui` is a Jupyter notebook GUI wrapper for the `ARMCore`
# simulator. It manages entering MOSFET parameters, and displaying
# relevent results and plots. This is built using `ipywidgets` and
# matplotlib ipywidgets backend package `ipympl`.
# 
# Original MATLAB ARM Core simulator written by Gage Hills and
# Prof. Max Shulaker for NOVELS group (2017).
# 
# Converted to Python by Andrew Yu <acyu@mit.edu> (NOVELS group)
# for 6.2500 ~4/20/2023.
 
import h5py
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

# CONSTANTS
LOG2 = 0.6931471805599453 # natural log(2)

def import_hdf5(path):
    """Import datasets in an hdf5 file into a dict.
    """
    d = {}

    with h5py.File(path, "r") as h5:
        for k in h5.keys():
            if len(h5[k].shape) > 0:
                d[k] = h5[k][:]
            else:
                d[k] = np.asscalar(h5[k][()])

    return d

def import_hdf5_dataset_group(group):
    """Import a flat hdf5 group into a dict. This does not work if the
    group contains nested groups.
    """
    d = {}

    for k in group.keys():
        if len(group[k].shape) > 0:
            d[k] = group[k][:]
        else:
            d[k] = np.asscalar(group[k][()])

    return d


def export_hdf5(path, data):
    """Export all keys in dict to hdf5 datasets, and save hdf5 file.
    """
    with h5py.File(path, "w") as h5:
        for k, val in data.items():
            h5.create_dataset(k, data=val)

@dataclass
class CellEnergyDelay:
    """Contains energy and delay parameters for a standard cell type,
    e.g. INV, NAND, etc.
    For notation of properties:
    - _r and _f: indicates rising and falling edge
    - _F, _A, _J etc.: indicates units of the value
    """
    name: str                               # cell name
    w_m: float                              # FET width [m]
    c_gs_F_m: float                         # gate to source capacitance per width [F/m]
    c_j_F_m: float                          # junction depletion capacitance at contacts per width [F/m]
    v_dd_V: float                           # supply voltage [V]
    i_on_A_m: float                         # on current [A/m]
    p_leak_baseline_nW: float               # leakage power at baseline [nW]
    num_fet_input: float                    # number of FETs on input, for calculating input cap
    num_fet_output: float                   # number of FETs on output, for calculating output parasitic cap
    num_fet_r: float                        # number of FETs driving current on rising edge
    num_fet_f: float                        # number of FETs driving current on falling edge
    
    # calculated, or set by user later
    c_in_F: float  = field(init=False)      # input capacitance [F]
    c_par_F: float = field(init=False)      # parasitic capacitance at output [F]
    i_dr_r_A: float = field(init=False)     # drive current during rising edge [A]
    i_dr_f_A: float = field(init=False)     # drive current during falling edge [A]
    t_d_par_r_s: float = field(init=False)  # delay of parasitic cap during rising edge [s]
    t_d_par_f_s: float = field(init=False)  # delay of parasitic cap during falling edge [s]
    energy_par_J: float = field(init=False) # energy of parasitic cap during one clock cycle [J]
    energy_int_J: float = 0                 # internal energy of cell during one clock cycle [J]

    def __post_init__(self):
        """Calculate derived parameters after initialization."""
        # num fets on input * gate capacitance per FET
        self.c_in_F = self.num_fet_input * self.c_gs_F_m * self.w_m
        
        # num fets on output * parasitic capacitance per FET
        self.c_par_F = self.num_fet_output * self.c_j_F_m * self.w_m

        # cell drive current during rising and falling edge [A]
        self.i_dr_r_A = self.num_fet_r * self.i_on_A_m * self.w_m # FETs driving pull-up
        self.i_dr_f_A = self.num_fet_f * self.i_on_A_m * self.w_m # FETs driving pull-down

        # cell output parasitic delay, from td ~ log(2) * Cpar * VDD / Ion
        self.t_d_par_r_s = LOG2 * self.c_par_F * self.v_dd_V / self.i_dr_r_A
        self.t_d_par_f_s = LOG2 * self.c_par_F * self.v_dd_V / self.i_dr_f_A

        # cell parasitic dynamic energy = 1/2 * Cpar * VDD^2
        self.energy_par_J = 0.5 * self.c_par_F * (self.v_dd_V**2)


class ARMCore():
    """Wrapper class for loading ARM core cells layout data.
    (So we only need to load data once.)
    With ARM core class instantiated, we calculate its energy and delay
    using `calculate_edp` method, which takes in configurable MOSFET
    parameters. General usage is something like:
    ```
    arm_core = ARMCore()
    metrics = arm_core.calculate_edp(
        i_on=1200,
        i_off=1,
        v_dd=1.8,
        cgs=1.2,
    )
    ```
    If inputs `i_on`, `i_off`, etc. are numpy vectors, this form allows
    quickly sweeping and calculating energy delay for a range of values.
    """

    def __init__(
        self,
        path_layout="armcore_layout.h5"
    ):
        """Loads ARM core cells from processed DEF and LEF file layout data,
        then converts relevent data into format easily usable by numpy.
        Instantiates a re-usable class container for the ARM core data and
        helpers for calculating energy and delay.
        """

        # baseline ARM core reference values
        # used for calculating leakage relative to baseline
        # DO NOT EDIT THESE
        self.baseline_v_dd_V = 0.5        # supply voltage [V]
        self.baseline_i_off_nA_um = 100.0 # off current [nA/um]
        self.baseline_w_nm = 72.0         # fet width [nm]
        self.baseline_w_m = 72e-9         # fet width [m]
        self.baseline_cpp_nm = 42.0       # contacted-poly-pitch [nm]

        # baseline gate and junction capacitance per FET
        # rough approximation for a high-k metal gate
        # DO NOT EDIT THESE
        self.EPS0   = 8.85e-12 # vacuum permittivity eps0 [F/m]
        self.H_GATE = 64e-9    # gate height [m]
        self.L_SPA  = 32e-9    # spacer length [m] 
        self.L_GATE = 32e-9    # gate length [m]
        self.L_C    = 32e-9    # contact length [m]
        self.XJ     = 4e-9     # junction depletion width [m]
        self.T_OX   = 2e-9     # oxide thickness [m]
        # self.K_OX   = 10.3     # oxide dielectric, realistic hi-K
        # self.K_SPA  = 5.5      # spacer dielectric, realistic SiN spacer
        self.K_OX   = 3.9      # oxide dielectric, 6.2500 project SiO2 gate oxide
        self.K_SPA  = 3.9      # spacer dielectric, 6.2500 project SiO2 spacer
        self.K_SI   = 11.7     # silicon dielectric
        self.CPP = self.L_GATE + self.L_C + 2*self.L_SPA # contacted poly pitch [m]

        # gate to channel capacitance per width [F/m]
        self.c_gc_F_m = self.L_GATE * self.K_OX/self.T_OX * self.EPS0
        # gate to source parasitic spacer capacitance per width [F/m]
        self.c_gs_par_F_m = self.H_GATE * self.K_SPA/self.L_SPA * self.EPS0
        # TOTAL gate to source capacitance per width [F/m]
        self.c_gs_F_m = self.c_gc_F_m + self.c_gs_par_F_m

        # junction depletion capacitance at contacts per width [F/m]
        self.c_j_F_m = self.L_C * self.K_SI/self.XJ * self.EPS0

        # LOAD ARM CORE LAYOUT DATA FROM h5 FILE
        with h5py.File(path_layout, "r") as h5:
            # parse attributes
            self.design_name =  h5.attrs["design_name"]
            self.diearea_xmin = h5.attrs["diearea_xmin"]
            self.diearea_ymin = h5.attrs["diearea_ymin"]
            self.diearea_xmax = h5.attrs["diearea_xmax"]
            self.diearea_ymax = h5.attrs["diearea_ymax"]
            self.units =        h5.attrs["units"]
            self.units_scale =  h5.attrs["units_scale"]

            # parse main data
            self.macros = import_hdf5_dataset_group(h5["macros"])
            self.components = import_hdf5_dataset_group(h5["components"])
            self.nets = import_hdf5_dataset_group(h5["nets"])
        
        ### DEBUG
        # print("self.macros", self.macros)
        # print(self.components["lib_cell"][0:10])
        # print(self.components["lib_cell_index"][0:10])
        # print(self.components["lib_cell_index"].dtype)
        
        # map a macro (standard lib cell) name to its macros index
        self.macro_name_to_index = {}
        for i, name_bytes in enumerate(self.macros["name"]):
            self.macro_name_to_index[name_bytes.decode("utf-8")] = i
        
        # lazy load polygons when layout plot is needed
        self.polygons = None
        self.polygons_dff = None
        self.polygons_inv_buf = None
        self.polygons_comb = None
    
    def generate_polygons(self):
        """Generate matplotlib polygon collections from layout data.
        """
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle
        
        idx_dff = self.components["index_DFF"]
        idx_inv_buf = np.logical_or(self.components["index_INV"], self.components["index_BUF"])
        idx_comb = np.logical_or(self.components["index_ND2"], self.components["index_NR2"])

        polygons = []
        for rect in self.components["layout_rect"]:
            polygons.append(Rectangle(
                (rect[0], rect[1]), # origin (x0, y0)
                rect[2], # width
                rect[3], # height
            ))
        self.polygons = np.array(polygons) # to allow fancy indexing

        self.polygons_dff = PatchCollection(
            self.polygons[idx_dff],
            facecolors='none',
            edgecolors=[0.85, 0, 0],
            antialiased=False,
            linewidths=0.5,
        )
        self.polygons_inv_buf = PatchCollection(
            self.polygons[idx_inv_buf],
            facecolors='none',
            edgecolors=[0, 0, 1],
            antialiased=False,
            linewidths=0.5,
        )
        self.polygons_comb = PatchCollection(
            self.polygons[idx_comb],
            facecolors='none',
            edgecolors=[0, 0.85, 0],
            antialiased=False,
            linewidths=0.5,
        )

    def create_cell(
        self,
        name: str,
        i_on_A_m: float,
        v_dd_V: float,
        c_gs: float, # normalized to 1
        w_m: float,
        p_leak_baseline_nW: float,
        num_fet_input: float,
        num_fet_output: float,
        num_fet_r: float,
        num_fet_f: float,
    ):
        """Convenience helper to create a standard cell energy delay
        CellEnergyDelay that automatically fills shared properties across
        arm core (w_m, c_gs, c_j).
        """
        return CellEnergyDelay(
            name=name,
            w_m=w_m,
            c_gs_F_m=c_gs * self.c_gs_F_m, # scale both c_gs and c_j by normalized
            c_j_F_m=c_gs * self.c_j_F_m,   # c_gs parameter
            v_dd_V=v_dd_V,
            i_on_A_m=i_on_A_m,
            p_leak_baseline_nW=p_leak_baseline_nW,
            num_fet_input=num_fet_input,
            num_fet_output=num_fet_output,
            num_fet_r=num_fet_r,
            num_fet_f=num_fet_f,
        )

    def calculate_energy_delay(
        self,
        
        # MOSFET configurable parameters
        i_on=600,           # in uA/um
        i_off=400,          # in nA/um
        v_dd=1.8,           # supply voltage
        c_gs=1,             # normalized to 1

        # built-in ARM core parameters
        # DO NOT EDIT THESE

        w_wire_m=32e-9,              # wire width [m]
        ar_wire=2,                   # wire aspect ratio
        c_wire_F_per_m = 0.15e-15 * 1e6, # wire cap per length in [F/m], default 0.15 fF/um
        rho_wire_Ohm_per_m=1.68e-8,      # wire resistivity, assuming copper [ohm*m]

        # FET width 
        # 7->10->14->22->32
        # scale by sqrt(2)^5 = 5.6529
        # 72 nm * 5.629 = 407.2935
        # 64 * ceil(72 * 5.629 / 64)
        w_nm = 448,          # fet width [nm]
        cpp_nm = 42,         # contacted poly pitch [nm]
        activity_factor=0.1, # dynamic power activity factor
        logic_depth=179,     # logic depth for ARM core
    ):
        """
        Returns
        results dict containing:
        - delay: float
        - energy_dyn: float
        - energy_static: float

        The general way this works:
        1. We have a list of standard cells and their energy delay parameters.
        2. We have list of components and nets in the design, which reference
           the standard cells.
        3. The nets form a graph between components. Unfortunately, this
           requires variable length arrays of drive components in each net.
        
        The storage format for current, capacitance, energy, delay, etc.
        parameters for calculations are struct-of-arrays style whenever
        possible, to maximize time in numpy vectorized operations and minimize
        python side loops.

        Macros: contains cap, energy, delay, etc. parameters
                 INDEX       0         1         2         3
                  name = [ INVD0 ] [ INVD1 ] [ ND2D0 ] [ BUFFD1 ] ...
                p_leak = [ 2.4   ] [ 3.6   ] [ 2.7   ] [ 7.2    ] ...
                c_in_F = [ 0.1   ] [ 0.2   ] [ 0.4   ] [ 0.8    ] ...
            delay_rise = [ 0.1   ] [ 0.1   ] [ 0.1   ] [ 0.2    ] ...
            delay_fall = [ 0.1   ] [ 0.1   ] [ 0.1   ] [ 0.2    ] ...

        Components/nets: index pointers reference index of library cells
        in `macros`:
            components = {
                # index into macros:
                lib_cell = [ 0 ] [ 0 ] [ 1 ] [ 1 ] [ 2 ] [ 3 ] ...
            }

        We can then use scatter/gather array indexing to get the parameters.
        E.g. for leakage power calculation, for each component we can do:

            components_p_leak = macros_p_leak[components_lib_cell]
            p_leak_total = np.sum(components_p_leak)
        
        For each component lib cell index, we index into the macros array to get
        its leakage. We can then sum up the leakage power for all components.

        - Converted from original Gage Hills / Max Shulaker MATLAB code by
        Andrew Yu for 6.2500 2023. 
        """
        # max int32 value, used as a sentinel for invalid index pointer
        MAX_INT = np.iinfo(np.int32).max

        # fet width in meters
        w_m = w_nm * 1e-9

        # wire derived parameters
        h_wire_m = ar_wire * w_wire_m  # wire height
        r_wire_Ohm_per_m = rho_wire_Ohm_per_m / (w_wire_m * h_wire_m) # wire resistance per length [Ohm/m]

        # standard cell energy delay parameters using new FET parameters
        cell_INVD0 = self.create_cell(
            name="INVD0",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=2.43185,
            num_fet_input=2*(2/3),  # 2 * (2/3) FETs on input
            num_fet_output=2*(2/3), # 2 * (2/3) FETs on output
            num_fet_r=2/3,          # 1 * (2/3) FETs full-drive pull-up
            num_fet_f=2/3,          # 1 * (2/3) FETs half-drive pull-down
        )

        cell_INVD1 = self.create_cell(
            name="INVD1",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=3.60024,
            num_fet_input=2,  # 2 FETs on input
            num_fet_output=2, # 2 FETs on output
            num_fet_r=1,      # 1 FET full-drive pull-up
            num_fet_f=1,      # 1 FET full-drive pull-down
        )

        cell_INVD2 = self.create_cell(
            name="INVD2",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=7.20049,
            num_fet_input=4,  # 4 FETs on input
            num_fet_output=4, # 4 FETs on output
            num_fet_r=2,      # 2 FETs full-drive pull-up
            num_fet_f=2,      # 2 FETs full-drive pull-down
        )

        cell_INVD4 = self.create_cell(
            name="INVD4",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=14.4009,
            num_fet_input=8,  # 8 FETs on input
            num_fet_output=8, # 8 FETs on output
            num_fet_r=4,      # 4 FETs full-drive pull-up
            num_fet_f=4,      # 4 FETs full-drive pull-down
        )

        cell_ND2D0 = self.create_cell(
            name="ND2D0",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=2.71848,
            num_fet_input=2 * (2/3),   # 2 * (2/3) FETs on input
            num_fet_output=3 * (2/3),  # 3 * (2/3) FETs on output
            num_fet_r=1 * (2/3),       # 1 * (2/3) FET full drive on pull-up
            num_fet_f=1 * (2/3) * 0.5, # 1 * (2/3) FET half drive on pull-down
        )

        cell_ND2D1 = self.create_cell(
            name="ND2D1",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=4.02458,
            num_fet_input=2,   # 2 FETs on input
            num_fet_output=3,  # 3 FETs on output
            num_fet_r=1,       # 1 FET full drive on pull-up
            num_fet_f=1 * 0.5, # 1 FET half drive on pull-down
        )

        cell_NR2ND2D0 = self.create_cell(
            name="NR2ND2D0",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=4.02458,
            num_fet_input=2 * (2/3),   # 2 * (2/3) FETs on input
            num_fet_output=3 * (2/3),  # 3 * (2/3) FETs on output
            num_fet_r=1 * (2/3) * 0.5, # 1 FET half drive on pull-up
            num_fet_f=1 * (2/3),       # 1 FET full drive on pull-down
        )

        # BUFFD1 is INVD1 driving INVD1
        cell_BUFFD1 = self.create_cell(
            name="BUFFD1",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=7.20137,
            num_fet_input=2,  # 2 FETs on input
            num_fet_output=2, # 2 FETs on output
            num_fet_r=1,      # 1 FET full-drive pull-up
            num_fet_f=1,      # 1 FET full-drive pull-down
        )
        cell_BUFFD1.t_d_par_r_s = (
            cell_INVD1.t_d_par_f_s                                        # INVD1 parasitic fall delay
            + LOG2 * cell_INVD1.c_in_F * v_dd / cell_INVD1.i_dr_f_A    # INVD1 fall to drive INVD1
            + LOG2 * cell_BUFFD1.c_par_F * v_dd / cell_BUFFD1.i_dr_r_A # BUFFD1 rise
        )
        cell_BUFFD1.t_d_par_f_s = (
            cell_INVD1.t_d_par_r_s                                        # INVD1 parasitic rise delay
            + LOG2 * cell_INVD1.c_in_F * v_dd / cell_INVD1.i_dr_r_A    # INVD1 rise to drive INVD1
            + LOG2 * cell_BUFFD1.c_par_F * v_dd / cell_BUFFD1.i_dr_f_A # BUFFD1 fall
        )
        cell_BUFFD1.energy_par_J = 0.5 * cell_BUFFD1.c_par_F * (v_dd**2)
        cell_BUFFD1.energy_int_J = (
            cell_INVD1.energy_par_J               # output of INVD1
            + 0.5 * cell_INVD1.c_in_F * (v_dd**2) # input of INVD1
        )

        # BUFFD2 is INVD1 driving INVD2
        cell_BUFFD2 = self.create_cell(
            name="BUFFD2",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=14.4028,
            num_fet_input=2,  # 2 FETs on input
            num_fet_output=4, # 4 FETs on output
            num_fet_r=2,      # 2 FET full-drive pull-up
            num_fet_f=2,      # 2 FET full-drive pull-down
        )
        cell_BUFFD2.t_d_par_r_s = (
            cell_INVD1.t_d_par_f_s                                     # INVD1 parasitic fall delay
            + LOG2 * cell_INVD2.c_in_F * v_dd / cell_INVD1.i_dr_f_A    # INVD1 fall to drive INVD2
            + LOG2 * cell_BUFFD2.c_par_F * v_dd / cell_BUFFD2.i_dr_r_A # BUFFD1 rise
        )
        cell_BUFFD2.t_d_par_f_s = (
            cell_INVD1.t_d_par_r_s                                     # INVD1 parasitic rise delay
            + LOG2 * cell_INVD1.c_in_F * v_dd / cell_INVD1.i_dr_r_A    # INVD1 rise to drive INVD1
            + LOG2 * cell_BUFFD2.c_par_F * v_dd / cell_BUFFD2.i_dr_f_A # BUFFD1 fall
        )
        cell_BUFFD2.energy_par_J = 0.5 * cell_BUFFD1.c_par_F * (v_dd**2)
        cell_BUFFD2.energy_int_J = (
            cell_INVD1.energy_par_J          # output of INVD1
            + 0.5 * cell_INVD2.c_in_F * (v_dd**2) # input of INVD2
        )

        # assume DFF has same input/output as BUFFD1
        # (approximate internal energy & setup/hold timing as 0)
        cell_DFCNQSTKND2D1_MOD3 = self.create_cell(
            name="DFCNQSTKND2D1_MOD3",
            i_on_A_m=i_on,
            v_dd_V=v_dd,
            c_gs=c_gs,
            w_m=w_m,
            p_leak_baseline_nW=26.7588,
            num_fet_input=cell_BUFFD2.num_fet_input,
            num_fet_output=cell_BUFFD2.num_fet_output,
            num_fet_r=cell_BUFFD2.num_fet_r,
            num_fet_f=cell_BUFFD2.num_fet_f,
        )
        cell_DFCNQSTKND2D1_MOD3.t_d_par_r_s = cell_BUFFD2.t_d_par_r_s
        cell_DFCNQSTKND2D1_MOD3.t_d_par_f_s = cell_BUFFD2.t_d_par_r_s # NOT typo, original code had this... (weird...)
        cell_DFCNQSTKND2D1_MOD3.energy_par_J = cell_BUFFD2.energy_par_J
        cell_DFCNQSTKND2D1_MOD3.energy_int_J = (
            9 * cell_ND2D0.energy_par_J +                  # 9 NAND outputs
            4 * cell_INVD0.energy_par_J +                  # 4 INV outputs
            0.5 * 9 * 2 * cell_ND2D0.c_in_F * (v_dd**2) +  # 9 NAND inputs (*2 inputs each)
            0.5 * 4 * cell_ND2D0.c_in_F * (v_dd**2)        # 4 INV inputs
        )

        # 1. create dict of cell.name => cell
        # 2. create arrays of relevent cell parameters
        #    (cap, delay, power, energy, etc.)
        #    where index corresponds to same cell in `self.macros`

        n_macros = self.macros["name"].shape[0]
        cells = {} # dict of cell.name => cell
        cell_p_leak_baseline_nW = np.zeros((n_macros,)) # leakage power at baseline [nW]
        cell_c_in_F  = np.zeros((n_macros,))            # input capacitance [F]
        cell_c_par_F = np.zeros((n_macros,))            # parasitic capacitance at output [F]
        cell_i_dr_r_A = np.zeros((n_macros,))           # drive current during rising edge [A]
        cell_i_dr_f_A = np.zeros((n_macros,))           # drive current during falling edge [A]
        cell_t_d_par_r_s = np.zeros((n_macros,))           # delay of parasitic cap during rising edge [s]
        cell_t_d_par_f_s = np.zeros((n_macros,))           # delay of parasitic cap during falling edge [s]
        cell_energy_par_J = np.zeros((n_macros,))       # energy of parasitic cap during one clock cycle [J]
        cell_energy_int_J = np.zeros((n_macros,))       # internal energy of cell during one clock cycle [J]

        for cell in [
            cell_BUFFD1,
            cell_BUFFD2,
            cell_DFCNQSTKND2D1_MOD3,
            cell_INVD0,
            cell_INVD1,
            cell_INVD2,
            cell_INVD4,
            cell_ND2D0,
            cell_ND2D1,
            cell_NR2ND2D0,
        ]:
            cells[cell.name] = cell

            idx = self.macro_name_to_index[cell.name]
            cell_p_leak_baseline_nW[idx] = cell.p_leak_baseline_nW
            cell_c_in_F[idx] = cell.c_in_F
            cell_c_par_F[idx] = cell.c_par_F
            cell_i_dr_r_A[idx] = cell.i_dr_r_A
            cell_i_dr_f_A[idx] = cell.i_dr_f_A
            cell_t_d_par_r_s[idx] = cell.t_d_par_r_s
            cell_t_d_par_f_s[idx] = cell.t_d_par_f_s
            cell_energy_par_J[idx] = cell.energy_par_J
            cell_energy_int_J[idx] = cell.energy_int_J

        # -----------------------------------------------------------
        # COMPUTE DELAY
        # -----------------------------------------------------------
        # accumulate load cap of each net that it is driving:
        # each net contains a variable length array of indices pointers into
        # cell macros indicating a list of load cells. we need to sum together
        # the input capacitances of all load cells to get total load capacitance
        # the net is driving.
        c_load_F = np.zeros((self.nets["name"].shape[0],), dtype=np.float64)
        for i, load_cell_indices in enumerate(self.nets["p_GATE_PAIR_lib_cell_load"]):
            c_load_F[i] = np.sum(cell_c_in_F[load_cell_indices])

        # contains net's input gate cell lib cell index 
        indices_drive_cells = self.nets["p_GATE_PAIR_lib_cell_drive_index"]
        valid_cells = indices_drive_cells != MAX_INT # invalid nets have no gates, e.g. some pin nets

        # cell drive current i_dr: use scatter/gather indexing into cell macro
        # i_dr arrays to gather i_dr_r and i_dr_f for each net
        i_dr_r_A = np.zeros_like(indices_drive_cells, dtype=np.float64)
        i_dr_r_A[valid_cells] = cell_i_dr_r_A[indices_drive_cells[valid_cells]]
        a = cell_i_dr_r_A[indices_drive_cells[valid_cells]]
        i_dr_r_A[valid_cells] = a
        i_dr_r_A[valid_cells][0] = a[0]
        i_dr_r_A[0] = a[0]

        i_dr_f_A = np.zeros_like(indices_drive_cells, dtype=np.float64)
        i_dr_f_A[valid_cells] = cell_i_dr_f_A[indices_drive_cells[valid_cells]]
        # wire resistance and capacitances load of each net
        # wire len:
        # scale wires in y dimension by the relative height of the standard cells
        wire_scale_y = w_nm / self.baseline_w_nm
        # scale wires in x dimension by the relative contacted gate pitch
        wire_scale_x = (self.CPP * 1e9) / self.baseline_cpp_nm
        wire_x_len_m = self.nets["x_len_m"] * wire_scale_x
        wire_y_len_m = self.nets["y_len_m"] * wire_scale_y
        wire_len_m = wire_x_len_m + wire_y_len_m # manhattan distance
        c_wire_F = wire_len_m * c_wire_F_per_m
        r_wire_Ohm = wire_len_m * r_wire_Ohm_per_m

        # cell parasitic delay: use scatter/gather indexing into cell macro
        # delay arrays to gather t_d_par delay for each net
        t_d_par_r_s = np.zeros_like(indices_drive_cells, dtype=np.float64)
        t_d_par_r_s[valid_cells] = cell_t_d_par_r_s[indices_drive_cells[valid_cells]]
        t_d_par_f_s = np.zeros_like(indices_drive_cells, dtype=np.float64)
        t_d_par_f_s[valid_cells] = cell_t_d_par_f_s[indices_drive_cells[valid_cells]]
        # wire elmore delays
        t_d_rwire_cwire_r_s = LOG2 * r_wire_Ohm * (c_wire_F/2)
        t_d_rwire_cwire_f_s = LOG2 * r_wire_Ohm * (c_wire_F/2)
        
        t_d_rdrive_cwire_r_s = LOG2 * np.divide( # on both sides of pi model
            (2 * c_wire_F / 2) * v_dd,
            i_dr_r_A,
            out=np.zeros_like(i_dr_r_A, dtype=np.float64),
            where=valid_cells, # avoids divide by zero
        )
        t_d_rdrive_cwire_f_s = LOG2 * np.divide( # on both sides of pi model
            (2 * c_wire_F / 2) * v_dd,
            i_dr_f_A,
            out=np.zeros_like(i_dr_f_A, dtype=np.float64),
            where=valid_cells, # avoids divide by zero
        )
        
        t_d_rwire_cload_r_s = LOG2 * r_wire_Ohm * c_load_F
        t_d_rwire_cload_f_s = LOG2 * r_wire_Ohm * c_load_F

        t_d_rdrive_cload_r_s = LOG2 * np.divide(
            c_load_F * v_dd,
            i_dr_r_A,
            out=np.zeros_like(i_dr_r_A, dtype=np.float64),
            where=valid_cells, # avoids divide by zero
        )
        t_d_rdrive_cload_f_s = LOG2 * np.divide(
            c_load_F * v_dd,
            i_dr_f_A,
            out=np.zeros_like(i_dr_f_A, dtype=np.float64),
            where=valid_cells, # avoids divide by zero
        )

        # rise and fall times of each component
        # = elmore delay of gate driving wire and fanout
        # (using pi model for distributed wire RC parasitics,
        # hence the rdrive*cwire and rwire*cin terms)
        t_d_r_all_s = (
            t_d_par_r_s
            + t_d_rwire_cwire_r_s
            + t_d_rdrive_cwire_r_s
            + t_d_rwire_cload_r_s
            + t_d_rdrive_cload_r_s
        )
        t_d_f_all_s = (
            t_d_par_f_s
            + t_d_rwire_cwire_f_s
            + t_d_rdrive_cwire_f_s
            + t_d_rwire_cload_f_s
            + t_d_rdrive_cload_f_s
        )

        ### DEBUG
        # print("t_d_par_r_s", t_d_par_r_s[:10])
        # print("t_d_par_f_s", t_d_par_f_s[:10])
        # print("t_d_rwire_cwire_r_s", t_d_rwire_cwire_r_s[:10])
        # print("t_d_rwire_cwire_f_s", t_d_rwire_cwire_f_s[:10])
        # print("t_d_rdrive_cwire_r_s", t_d_rdrive_cwire_r_s[:10])
        # print("t_d_rdrive_cwire_f_s", t_d_rdrive_cwire_f_s[:10])
        # print("t_d_rwire_cload_r_s", t_d_rwire_cload_r_s[:10])
        # print("t_d_rwire_cload_f_s", t_d_rwire_cload_f_s[:10])
        # print("t_d_rdrive_cload_r_s", t_d_rdrive_cload_r_s[:10])
        # print("t_d_rdrive_cload_f_s", t_d_rdrive_cload_f_s[:10])
        # print("t_d_r_all_s", t_d_r_all_s[:10])

        # filter only valid cells
        t_d_r_filtered_s = t_d_r_all_s[valid_cells]
        t_d_f_filtered_s = t_d_f_all_s[valid_cells]

        # average rising and falling edge delay of each net
        t_d_r_avg_s = np.mean(t_d_r_filtered_s)
        t_d_f_avg_s = np.mean(t_d_f_filtered_s)

        # take average of rising and falling edge delays
        t_d_avg_s = 0.5 *(t_d_r_avg_s + t_d_f_avg_s)

        # delay of critical path (of logic_depth number of stages)
        # = logic_depth * average delay of each stage
        t_crit_s = logic_depth * t_d_avg_s

        # clock frequency = 1/critical path delay
        clk_freq_Hz = 1 / t_crit_s
        clk_freq_GHz = clk_freq_Hz / 1e9

        # -----------------------------------------------------------
        # COMPUTE DYNAMIC ENERGY
        # -----------------------------------------------------------
        # total energy from parasitic capacitances in nets
        energy_par_J = np.zeros_like(indices_drive_cells, dtype=np.float64)
        energy_par_J[valid_cells] = cell_energy_par_J[indices_drive_cells[valid_cells]]

        # load capacitance energy per cell
        energy_load_J = 0.5 * c_load_F * (v_dd**2)
        
        # wire energy per cell
        energy_wire_J = 0.5 * c_wire_F * (v_dd**2)

        # total energy (in valid cells)
        energy_par_total_J = np.sum(energy_par_J[valid_cells])
        energy_load_total_J = np.sum(energy_load_J[valid_cells])
        energy_wire_total_J = np.sum(energy_wire_J)

        # dynamic energy per cycle ~ C * V^2 * f * activity_factory
        energy_cycle_par_total_J = activity_factor * energy_par_total_J
        energy_cycle_wire_total_J = activity_factor * energy_wire_total_J
        energy_cycle_load_total_J = activity_factor * energy_load_total_J
        energy_dyn_cycle_total_J = energy_cycle_par_total_J + energy_cycle_wire_total_J + energy_cycle_load_total_J

        # -----------------------------------------------------------
        # COMPUTE LEAKAGE POWER AND ENERGY PER CYCLE
        # -----------------------------------------------------------
        # maps each cell's macros array index in components to that cell's
        # leakage power parameter (using scatter/gather array indexing),
        # then sum together to get total leakage power
        p_leak_ref_total_nW = np.sum(cell_p_leak_baseline_nW[self.components["lib_cell_index"]])
        
        ### DEBUG
        # print(cell_p_leak_baseline_nW[self.components["lib_cell_index"]][0:20]) 

        # calculate leakage power by scaling baseline reference leakage
        # by new VDD parameter and Ioff parameter
        p_leak_total_W = (
            1e-9 * p_leak_ref_total_nW
            * (v_dd / self.baseline_v_dd_V) # scale by voltage
            * (i_off / self.baseline_i_off_nA_um) # scale by i_off
            * (w_nm / self.baseline_w_nm) # scale by width
        )

        # calculate leakage energy per cycle by multiplying leakage power
        # by cycle time (need to calculate cycle delay first)
        energy_leak_total_J = p_leak_total_W / clk_freq_Hz

        # -----------------------------------------------------------
        # COMPUTE TOTAL ENERGY PER CYCLE AND EDP
        # -----------------------------------------------------------
        # total energy per cycle
        energy_cycle_total_J = energy_dyn_cycle_total_J + energy_leak_total_J
        
        # energy-delay product
        edp_Js = energy_cycle_total_J * t_crit_s
        edp_pJns = edp_Js * 1e12 * 1e9

        # pre-convert energy to pJ units for easier display/plotting
        energy_cycle_par_total_pJ = 1e12 * energy_cycle_par_total_J
        energy_cycle_wire_total_pJ = 1e12 * energy_cycle_wire_total_J
        energy_cycle_load_total_pJ = 1e12 * energy_cycle_load_total_J
        energy_dyn_cycle_total_pJ = 1e12 * energy_dyn_cycle_total_J
        energy_leak_total_pJ = 1e12 * energy_leak_total_J
        energy_cycle_total_pJ = 1e12 * energy_cycle_total_J

        # return results as dict
        results = {
            "i_on_A_m": i_on,
            "i_off_nA_um": i_off,
            "v_dd_V": v_dd,
            "c_gs": c_gs,
            "t_crit_s": t_crit_s,
            "clk_freq_GHz": clk_freq_GHz,
            "p_leak_total_W": p_leak_total_W,
            # store energy in pJ for easier display/plotting
            "energy_cycle_par_total_pJ": energy_cycle_par_total_pJ,
            "energy_cycle_wire_total_pJ": energy_cycle_wire_total_pJ,
            "energy_cycle_load_total_pJ": energy_cycle_load_total_pJ,
            "energy_dyn_cycle_total_pJ": energy_dyn_cycle_total_pJ,
            "energy_leak_total_pJ": energy_leak_total_pJ,
            "energy_cycle_total_pJ": energy_cycle_total_pJ,
            "edp_Js": edp_Js,
            "edp_pJns": edp_pJns,
        }

        return results


class EDPGui():
    """Implements a Jupyter notebook ipywidgets GUI for inputting FET
    parameters and calculating EDP.
    """

    def __init__(
        self,
        path_results="results.csv",
        debug=False, # turn on debug print lines
    ):
        import os
        import traceback
        from IPython.display import display, HTML
        import ipywidgets as widgets

        # debug messages
        self.debug = debug

        # arm core simulation instance
        self.arm_core = ARMCore()

        # initial values in GUI
        self.initial_i_on = 600.0  # in uA/um
        self.initial_i_off = 400.0 # in nA/um
        self.initial_v_dd = 1.8  # supply voltage
        self.initial_c_gs = 1.0  # relative Cgs, normalized to 1
        
        # list of tracked result keys from simulation
        self.tracked_results = [
            "i_on_A_m",
            "i_off_nA_um",
            "v_dd_V",
            "c_gs",
            "t_crit_s",
            "clk_freq_GHz",
            "p_leak_total_W",
            "energy_cycle_par_total_pJ",
            "energy_cycle_wire_total_pJ",
            "energy_cycle_load_total_pJ",
            "energy_dyn_cycle_total_pJ",
            "energy_leak_total_pJ",
            "energy_cycle_total_pJ",
            "edp_Js",
            "edp_pJns",
        ]

        # string header for csv file from tracked result keys
        self.csv_tracked_results_header = ",".join(self.tracked_results)
        
        # store file path to store results
        self.path_results = path_results

        # store history vector for each tracked result key
        # inside this dict
        self.results_history = {}
        for key in self.tracked_results:
            self.results_history[key] = []

        # load stored results if not None
        loaded_stored_results = False
        if path_results is not None and os.path.exists(path_results):
            try:
                with open(path_results, "r") as f:
                    stored_results = f.readlines()
                if len(stored_results) > 1:
                    for line in stored_results[1:]: # skip header
                        line = line.strip().split(",")
                        for i, key in enumerate(self.tracked_results):
                            self.results_history[key].append(float(line[i]))
                    # set initial values to last stored values
                    self.initial_i_on = self.results_history["i_on_A_m"][-1]
                    self.initial_i_off = self.results_history["i_off_nA_um"][-1]
                    self.initial_v_dd = self.results_history["v_dd_V"][-1]
                    self.initial_c_gs = self.results_history["c_gs"][-1]
                    loaded_stored_results = True
            except Exception as e:
                print(repr(e))
                traceback.print_exc()

        # otherwise, do initial simulation and store initial results
        if not loaded_stored_results:
            results = self.arm_core.calculate_energy_delay(
                i_on=self.initial_i_on,
                i_off=self.initial_i_off,
                v_dd=self.initial_v_dd,
                c_gs=self.initial_c_gs,
            )
            for key in self.tracked_results:
                self.results_history[key].append(results[key])

        # =============================================================================
        # FET Parameters Inputs
        # =============================================================================
        self.label_fet_params = widgets.Label(
            value="MOSFET Parameters From Sentaurus",
            style=dict(
                font_weight="bold",
                font_size="16px",
            ),
            layout=widgets.Layout(
                font_size="48px",
                display="flex",
                flex="1 1 auto",
                justify_content="center",
                overflow="visible",
            ),
        )

        self.box_fet_param_title = widgets.HBox(
            children=[
                self.label_fet_params,
            ],
            layout=widgets.Layout(
                display="flex",
                flex="1 1 auto",
                flex_grow="1",
                overflow="visible",
            ),
        )

        # on current input
        self.input_fet_param_label_on_current = widgets.Label(
            value="On Current [µA/µm] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.input_fet_param_on_current = widgets.FloatText(
            value=self.initial_i_on,
            step=1,
            layout=widgets.Layout(
                width="35%",
            ),
        )

        self.box_fet_param_on_current = widgets.HBox(
            children=[
                self.input_fet_param_label_on_current,
                self.input_fet_param_on_current,
            ],
        )

        # off current input
        self.input_fet_param_label_off_current = widgets.Label(
            value="Off Current [nA/µm] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.input_fet_param_off_current = widgets.FloatText(
            value=self.initial_i_off,
            step=1,
            layout=widgets.Layout(
                width="35%",
            ),
        )

        self.box_fet_param_off_current = widgets.HBox(
            children=[
                self.input_fet_param_label_off_current,
                self.input_fet_param_off_current,
            ],
        )

        # supply voltage VDD input
        self.input_fet_param_label_vdd = widgets.Label(
            value="Supply Voltage [V] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.input_fet_param_v_dd = widgets.FloatText(
            value=self.initial_v_dd,
            step=0.001,
            layout=widgets.Layout(
                width="35%",
            ),
        )

        self.box_fet_param_vdd = widgets.HBox(
            children=[
                self.input_fet_param_label_vdd,
                self.input_fet_param_v_dd,
            ],
        )

        # gate-to-source capacitance Cgs input
        self.input_fet_param_label_cgs = widgets.VBox(
            children=[
                widgets.Label(
                    value="Gate-to-Source Capacitance = ",
                    style=dict(
                        font_size="13px",
                    ),
                    layout=widgets.Layout(
                        display="flex",
                        justify_content="flex-end",
                        padding="0",
                        margin="-4px 0 0 0",
                    ),
                ),
                widgets.Label(
                    value="(Normalized)",
                    style=dict(
                        font_size="13px",
                    ),
                    layout=widgets.Layout(
                        display="flex",
                        justify_content="flex-end",
                        padding="0 12px",
                        margin="-6px 0 0 0",
                    ),
                ),
            ],
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
                margin="2px",
            ),
        )

        self.input_fet_param_c_gs = widgets.FloatText(
            value=self.initial_c_gs,
            step=0.001,
            layout=widgets.Layout(
                width="35%",
            ),
        )

        self.box_fet_param_cgs = widgets.HBox(
            children=[
                self.input_fet_param_label_cgs,
                self.input_fet_param_c_gs,
            ],
        )

        # button to calculate edp with input values
        self.button_fet_param_calculate_edp = widgets.Button(
            description="Calculate EDP",
            layout=widgets.Layout(
                width="60%",
                border="1px solid #aaa",
            ),
        )

        self.box_fet_param_button_calculate_edp = widgets.HBox(
            children=[
                self.button_fet_param_calculate_edp,
            ],
            layout=widgets.Layout(
                width="100%",
                display="flex",
                justify_content="center",
            ),
        )
        
        # callback to calculate EDP when button pressed
        self.button_fet_param_calculate_edp.on_click(
            lambda _: self.calculate(),
        )

        # box containing fet parameters inputs
        self.box_fet_params = widgets.VBox(
            children=[
                self.box_fet_param_title,
                self.box_fet_param_on_current,
                self.box_fet_param_off_current,
                self.box_fet_param_vdd,
                self.box_fet_param_cgs,
                self.box_fet_param_button_calculate_edp,
            ],
            layout=widgets.Layout(
                border="1px solid #222222",
                margin="2px",
                padding="2px 0px 8px 0px",
            ),
        )

        # =============================================================================
        # Derived Metrics from FET Parameters
        # =============================================================================

        # first repeat the fet parameter input values
        # (to make sure the user does not accidentally change value
        # in inputs and not remember values used to calculate metrics)
        self.label_metrics_fet_params =  widgets.Label(
            value="MOSFET Parameters",
            style=dict(
                font_weight="bold",
                font_size="16px",
            ),
            layout=widgets.Layout(
                display="flex",
                justify_content="center",
            ),
        )

        # metric fet param: on current
        self.metrics_fet_on_current_label = widgets.Label(
            value="On Current [µA/µm] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
                background_color="#aaa",
            ),
        )

        self.metrics_fet_on_current = widgets.Label(
            value=f"{self.input_fet_param_on_current.value:.2f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_fet_on_current = widgets.HBox(
            children=[
                self.metrics_fet_on_current_label,
                self.metrics_fet_on_current,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric fet param: off current
        self.metrics_fet_off_current_label = widgets.Label(
            value="Off Current [nA/µm] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_fet_off_current = widgets.Label(
            value=f"{self.input_fet_param_off_current.value:.2f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_fet_off_current = widgets.HBox(
            children=[
                self.metrics_fet_off_current_label,
                self.metrics_fet_off_current,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric fet param: supply voltage VDD
        self.metrics_fet_v_dd_label = widgets.Label(
            value="Supply VDD [V] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_fet_v_dd = widgets.Label(
            value=f"{self.input_fet_param_v_dd.value:.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_fet_v_dd = widgets.HBox(
            children=[
                self.metrics_fet_v_dd_label,
                self.metrics_fet_v_dd,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric fet param: gate-to-source cap Cgs
        self.metrics_fet_c_gs_label = widgets.Label(
            value="Cgs (Normalized) = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_fet_c_gs = widgets.Label(
            value=f"{self.input_fet_param_c_gs.value:.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_fet_c_gs = widgets.HBox(
            children=[
                self.metrics_fet_c_gs_label,
                self.metrics_fet_c_gs,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # ARM CORE CALCULATED EDP METRICS
    
        self.label_metrics_calculated_edp =  widgets.Label(
            value="Calculated ARM Core Metrics",
            style=dict(
                font_weight="bold",
                font_size="16px",
            ),
            layout=widgets.Layout(
                display="flex",
                justify_content="center",
            ),
        )

        # metric: leakage energy
        self.metrics_energy_leakage_label = widgets.Label(
            value="Leakage Energy [pJ] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_energy_leakage = widgets.Label(
            value=f"{(self.results_history['energy_leak_total_pJ'][-1]):.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_energy_leakage = widgets.HBox(
            children=[
                self.metrics_energy_leakage_label,
                self.metrics_energy_leakage,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric: dynamic energy
        self.metrics_energy_dynamic_label = widgets.Label(
            value="Dynamic Energy [pJ] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_energy_dynamic = widgets.Label(
            value=f"{(self.results_history['energy_dyn_cycle_total_pJ'][-1]):.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_energy_dynamic = widgets.HBox(
            children=[
                self.metrics_energy_dynamic_label,
                self.metrics_energy_dynamic,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric: total energy
        self.metrics_energy_total_label = widgets.Label(
            value="TOTAL Energy [pJ] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_energy_total = widgets.Label(
            value=f"{(self.results_history['energy_cycle_total_pJ'][-1]):.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_energy_total = widgets.HBox(
            children=[
                self.metrics_energy_total_label,
                self.metrics_energy_total,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric: clock frequency
        self.metrics_freq_label = widgets.Label(
            value="Clock Frequency [GHz] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_freq = widgets.Label(
            value=f"{self.results_history['clk_freq_GHz'][-1]:.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_freq = widgets.HBox(
            children=[
                self.metrics_freq_label,
                self.metrics_freq,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # metric: edp
        self.metrics_edp_label = widgets.Label(
            value="EDP [pJ⋅ns] = ",
            layout=widgets.Layout(
                width="60%",
                display="flex",
                justify_content="flex-end",
            ),
        )

        self.metrics_edp = widgets.Label(
            value=f"{self.results_history['edp_pJns'][-1]:.3f}",
            layout=widgets.Layout(
                width="30%",
            ),
        )

        self.box_metrics_edp = widgets.HBox(
            children=[
                self.metrics_edp_label,
                self.metrics_edp,
            ],
            layout=widgets.Layout(
                overflow="hidden",
            ),
        )

        # box containing armcore metrics
        self.box_metrics = widgets.VBox(
            children=[
                self.label_metrics_fet_params,
                self.box_metrics_fet_on_current,
                self.box_metrics_fet_off_current,
                self.box_metrics_fet_v_dd,
                self.box_metrics_fet_c_gs,
                self.label_metrics_calculated_edp,
                self.box_metrics_energy_leakage,
                self.box_metrics_energy_dynamic,
                self.box_metrics_energy_total,
                self.box_metrics_freq,
                self.box_metrics_edp,
            ],
            layout=widgets.Layout(
                border="1px solid #222222",
                margin="2px",
            ),
        )

        # =============================================================================
        # Layout plot
        # =============================================================================
        self.label_layout =  widgets.Label(
            value="ARM Core Physical Layout",
            style=dict(
                font_weight="bold",
                font_size="16px",
            ),
            layout=widgets.Layout(
                display="flex",
                justify_content="center",
            ),
        )

        self.checkbox_layout_plot_inv_buf = widgets.Checkbox(
            value=True,
            indent=False,
            description="Show Inverters & Buffers",
            layout=widgets.Layout(
                padding="0px 0px 0px 40px",
            ),
        )
        self.checkbox_layout_plot_comb = widgets.Checkbox(
            value=True,
            indent=False,
            description="Show Combinational Logic",
            layout=widgets.Layout(
                padding="0px 0px 0px 40px",
            ),
        )
        self.checkbox_layout_plot_dff = widgets.Checkbox(
            value=True,
            indent=False,
            description="Show Sequential Logic: D-Flip-Flops",
            layout=widgets.Layout(
                padding="0px 0px 0px 40px",
            ),
        )

        self.checkbox_layout_plot_inv_buf.add_class("checkbox-label-blue")
        self.checkbox_layout_plot_comb.add_class("checkbox-label-green")
        self.checkbox_layout_plot_dff.add_class("checkbox-label-red")
        
        display(HTML("""<style>
        .checkbox-label-blue > .widget-label-basic {color: #0000DD; font-weight: bold;}
        .checkbox-label-green > .widget-label-basic {color: #00BB00; font-weight: bold;}
        .checkbox-label-red > .widget-label-basic {color: #DD0000; font-weight: bold;}
        </style>"""))

        self.output_plot_layout = widgets.Output(
            layout=widgets.Layout(
                overflow="visible",
            )
        )

        # create figure for physical core layout
        with self.output_plot_layout:
            self.fig_layout, self.ax_layout = plt.subplots(
                figsize=(3, 3),
            )
            self.fig_layout.canvas.header_visible = False # get rid of "Figure 1" title
            self.ax_layout.set_position([0.2, 0.2, 0.6, 0.6])
            plt.show()
        
        # do initial layout plot
        self.plot_layout()

        # handlers to replot layout when checkboxes are toggled
        self.checkbox_layout_plot_inv_buf.observe(
            lambda _: self.plot_layout(),
            names=["value"],
        )
        self.checkbox_layout_plot_comb.observe(
            lambda _: self.plot_layout(),
            names=["value"],
        )
        self.checkbox_layout_plot_dff.observe(
            lambda _: self.plot_layout(),
            names=["value"],
        )

        # box containing armcore metrics
        self.box_layout_plot = widgets.VBox(
            children=[
                self.label_layout,
                self.checkbox_layout_plot_inv_buf,
                self.checkbox_layout_plot_comb,
                self.checkbox_layout_plot_dff,
                self.output_plot_layout,
            ],
            layout=widgets.Layout(
                border="1px solid #222222",
                margin="2px",
            ),
        )

        # =============================================================================
        # EDP and energy breakdown plots
        # =============================================================================
        # plot for energy vs. frequency
        self.output_plot_energy_freq = widgets.Output(
            layout=widgets.Layout(),
        )
        with self.output_plot_energy_freq:
            self.fig_energy_freq, self.ax_energy_freq = plt.subplots(
                figsize=(5, 4),
            )
            self.fig_energy_freq.canvas.header_visible = False # get rid of "Figure 2" title
            plt.show()
        
        # plot for energy breakdown
        self.output_plot_energy_breakdown = widgets.Output(
            layout=widgets.Layout(),
        )

        with self.output_plot_energy_breakdown:
            self.fig_energy_breakdown, self.ax_energy_breakdown = plt.subplots(
                figsize=(5, 2),
            )
            self.ax_energy_breakdown.set_position([0.12, 0.3, 0.58, 0.5])
            self.fig_energy_breakdown.canvas.header_visible = False # get rid of "Figure 2" title
            plt.show()
        
        # do initial plot
        self.plot_energy_freq()

        # right side panel box containing output plots 
        self.box_right_plots = widgets.VBox(
            children=[
                self.output_plot_energy_freq, 
                self.output_plot_energy_breakdown,
            ],
            layout=widgets.Layout(
                border="1px solid #222222",
                margin="2px",
                overflow="scroll",
                flex="1 0 auto", # prevent shrink
            ),
        )

        # =============================================================================
        # Overall layout
        # =============================================================================
        # box containing fet parameters, fet metrics, and armcore metrics
        self.box_left_params = widgets.VBox(
            children=[self.box_fet_params, self.box_metrics, self.box_layout_plot],
            layout=widgets.Layout(
                width="auto",
                flex="1 0 auto", # prevent shrink
            )
        )

        # final box for app layout
        self.box_app = widgets.HBox(
            children=[self.box_left_params, self.box_right_plots],
            layout=widgets.Layout(
                width="auto",
                align_items="stretch",
                overflow="scroll",
            )
        )

    def plot_energy_freq(self):
        """Re-plot energy vs. clock frequency and energy breakdowns."""
        len_history = len(self.results_history["clk_freq_GHz"])

        with self.output_plot_energy_freq:
            self.ax_energy_freq.clear()

            self.ax_energy_freq.set_title(f"ARM Core: Energy vs. Clock Freq")
            self.ax_energy_freq.set_xlabel("Clock Frequency [GHz]")
            self.ax_energy_freq.set_ylabel("Energy [pJ]")
            self.ax_energy_freq.grid(
                color=[0.85, 0.85, 0.85],
            )

            h0 = self.ax_energy_freq.plot(
                self.results_history["clk_freq_GHz"],
                self.results_history["energy_cycle_total_pJ"],
                linewidth=1.0,
                color=[0, 0, 0],
                marker="o",
                markersize=4,
            )

            # plot single scatter on previous and most recent design
            if len_history >= 2: # previous design
                h1 = self.ax_energy_freq.scatter(
                    self.results_history["clk_freq_GHz"][-2],
                    self.results_history["energy_cycle_total_pJ"][-2],
                    color=[[0, 0, 1]],
                    marker="o",
                    s=60,
                    zorder=3,
                )
            
            if len_history >= 1: # current design
                h2 = self.ax_energy_freq.scatter(
                    self.results_history["clk_freq_GHz"][-1],
                    self.results_history["energy_cycle_total_pJ"][-1],
                    color=[[1, 0, 0]],
                    marker="o",
                    s=60,
                    zorder=3,
                )

            # legend
            if len_history >= 1:
                if len_history >= 2:
                    handles = (h0[0], h1, h2)
                    labels = ("All Designs", "Previous", "Current")
                elif len_history == 1:
                    handles = (h0[0], h2)
                    labels = ("All Designs", "Current")
                self.ax_energy_freq.legend(
                    handles=handles,
                    labels=labels,
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                )
            else:
                # no designs plotted
                pass

            self.fig_energy_freq.tight_layout()
        
        # energy breakdown of previous and current design
        with self.output_plot_energy_breakdown:
            self.ax_energy_breakdown.clear()
            
            # horizontal bar plot showing current and previous 
            # energy breakdown
            energy_and_color = [ # (energy_name, energy_key, colors)
                ("FET Cgate", "energy_cycle_load_total_pJ", "#1f77e4"),
                ("FET Cpar", "energy_cycle_par_total_pJ", "#2ca02c"),
                ("Wires", "energy_cycle_wire_total_pJ", "#d62728"),
                ("Leakage", "energy_leak_total_pJ", "#ff7f0e"),
            ]

            energy_curr_sum = 0
            energy_prev_sum = 0
            
            handles = []

            for (energy_name, energy_key, color) in energy_and_color:
                if len_history >= 1:
                    energy_curr_component = self.results_history[energy_key][-1]
                    
                    h = self.ax_energy_breakdown.barh(
                        "Curr",
                        [energy_curr_component],
                        height=0.75,
                        label=energy_name,
                        left=energy_curr_sum,
                        color=color,
                    )
                    handles.append(h)

                    energy_curr_sum += energy_curr_component

                if len_history >= 2:
                    energy_prev_component = self.results_history[energy_key][-2]

                    self.ax_energy_breakdown.barh(
                        "Prev",
                        [energy_prev_component],
                        height=0.75,
                        label=energy_name,
                        left=energy_prev_sum,
                        color=color,
                    )

                    energy_prev_sum += energy_prev_component
                
            # only adds legend for one bar plot to reduce clutter
            # (color is same for both plots)
            if len(handles) > 0:
                self.ax_energy_breakdown.legend(
                    handles=handles,
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                )
                
            self.ax_energy_breakdown.set_title(f"ARM Core: Energy Breakdown")
            self.ax_energy_breakdown.set_xlabel("Energy [pJ]")

    def plot_layout(self):
        """Plot arm core physical layout."""
        if self.arm_core.polygons is None:
            self.arm_core.generate_polygons()
        
        with self.output_plot_layout:
            self.ax_layout.clear()

            self.ax_layout.set_title(f"ARM Core Physical Layout")
            self.ax_layout.set_xlabel("x [µm]")
            self.ax_layout.set_ylabel("y [µm]")
            self.ax_layout.set_xlim(self.arm_core.diearea_xmin - 1, self.arm_core.diearea_xmax + 1)
            self.ax_layout.set_ylim(self.arm_core.diearea_ymin - 1, self.arm_core.diearea_ymax + 1)
            self.ax_layout.set_aspect("equal", "box")
            
            if self.checkbox_layout_plot_comb.value == True:
                self.ax_layout.add_collection(self.arm_core.polygons_comb)
            if self.checkbox_layout_plot_inv_buf.value == True:
                self.ax_layout.add_collection(self.arm_core.polygons_inv_buf)
            if self.checkbox_layout_plot_dff.value == True:
                self.ax_layout.add_collection(self.arm_core.polygons_dff)
    
    def save_results_to_csv(self):
        # number of history items (get length of any tracked result vector)
        n_history = len(self.results_history[self.tracked_results[0]])
        
        # number of tracked items
        n_items = len(self.tracked_results)

        # allocate list of `n_items` sized empty lists representing each history line
        history_lines = [[None] * n_items for i in range(n_history)]
        
        # populate history lines
        for i, key in enumerate(self.tracked_results):
            for j, val in enumerate(self.results_history[key]):
                history_lines[j][i] = str(val)
        
        with open(self.path_results, "w+") as f:
            f.write(self.csv_tracked_results_header + "\n")
            # write history lines
            for vals in history_lines:
                f.write(",".join(vals) + "\n")

    def calculate(self):
        """Run EDP calculation:
        1. Run EDP calculation
        2. Update metrics display on right side panel
        3. Update EDP frequency vs energy plot
        """
        # gather MOSFET input parameters
        i_on = self.input_fet_param_on_current.value
        i_off = self.input_fet_param_off_current.value
        v_dd = self.input_fet_param_v_dd.value
        c_gs = self.input_fet_param_c_gs.value

        if self.debug:
            print("Calculating EDP...")
            print(f"i_on = {i_on}, i_off = {i_off}, vdd = {v_dd}, cgs = {c_gs}")
        
        # run EDP calculation
        results = self.arm_core.calculate_energy_delay(
            i_on=i_on,
            i_off=i_off,
            v_dd=v_dd,
            c_gs=c_gs,
        )

        if self.debug:
            print(f"results: {results}")

        # update metrics display on right side panel
        self.metrics_fet_on_current.value = f"{i_on:.2f}"
        self.metrics_fet_off_current.value = f"{i_off:.2f}"
        self.metrics_fet_v_dd.value = f"{v_dd:.3f}"
        self.metrics_fet_c_gs.value = f"{c_gs:.3f}"
        self.metrics_energy_leakage.value = f"{(results['energy_leak_total_pJ']):.3f}"
        self.metrics_energy_dynamic.value = f"{(results['energy_dyn_cycle_total_pJ']):.3f}"
        self.metrics_energy_total.value = f"{(results['energy_cycle_total_pJ']):.3f}"
        self.metrics_freq.value = f"{results['clk_freq_GHz']:.3f}"
        self.metrics_edp.value = f"{results['edp_pJns']:.3f}"

        # store results into history and save history
        for key in self.tracked_results:
            self.results_history[key].append(results[key])
        if self.path_results is not None:
            self.save_results_to_csv()
        
        # re-plot
        self.plot_energy_freq()
    
    def display(self):
        """Displays gui in notebook."""
        from IPython.display import display, clear_output
        display(self.box_app)


if __name__ == "__main__":
    # FOR NOW: command line operation will just do verification run
    # to check if values match old MATLAB simulation values

    # create an instance of the ARMCore simulator
    armcore = ARMCore()

    print("Running ARM Core EDP Verification with MATLAB results...")

    # verify a few values
    for (i_on, i_off, v_dd, c_gs, total_energy_expected, clock_freq_expected, edp_expected) in [
        (600, 400, 1.8, 1.0,  32.14,  0.2401, 133.8),
        (1000, 400, 1.8, 1.0, 22.45,  0.3954,  56.78),
        (1000, 10, 1.8, 1.0,   7.854, 0.3954,  19.86),
        (600, 400, 1.0, 1.0,  10.03,  0.4260,  23.54),
        (600, 400, 1.8, 0.5,  23.47,  0.3317,  70.77),
        (1000, 100, 1.0, 0.5,  2.595, 0.9556,   2.716),
        (1000, 100, 3.0, 2.0, 47.04,  0.1547, 304.13),
    ]:
        # calculate edp
        results = armcore.calculate_energy_delay(
            i_on=i_on,     # in uA/um
            i_off=i_off,   # in nA/um
            v_dd=v_dd,     # supply voltage
            c_gs=c_gs,     # normalized to 1
        )

        # difference in %
        diff_total_energy = 100.0 * (total_energy_expected - results["energy_cycle_total_pJ"])
        diff_clk_freq = 100.0 * (clock_freq_expected - results["clk_freq_GHz"])
        diff_edp = 100.0 * (edp_expected - results["edp_pJns"])

        print("================================================================")
        print(f"i_on = {i_on}, i_off = {i_off}, vdd = {v_dd}, cgs = {c_gs}")
        print("PARAMETER: EXPECTED | CALCULATED | DIFFERENCE [%]")
        print(f"total_energy [pJ]: {total_energy_expected} | {results['energy_cycle_total_pJ']} | {diff_total_energy:.2f}")
        print(f"clock_freq [GHz]: {clock_freq_expected} | {results['clk_freq_GHz']} | {diff_clk_freq:.2f}")
        print(f"edp [pJ*ns]: {edp_expected} | {results['edp_pJns']} | {diff_edp:.2f}")
