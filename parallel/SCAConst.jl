module SCAConst

export PenConflict, PenMinSep, PenCloseness, PenAction, StdDist, MinSepSq, InvVar
export DT, DTI, G, Xmin, Xmax, Ymin, Ymax, Pmin, Pmax, Vmin, Vmax
export Xdim, Ydim, Pdim, Vdim, COCdim, NStates, NActions

const PenConflict = 1.0
const PenMinSep = 1000.0
const PenCloseness = 10.0
const PenAction = 0.02

const StdDist = 500.0  # [m]
const MinSepSq = StdDist^2  # [m^2]
const InvVar = 1 / MinSepSq  # [1/m^2]

const DT = 5.0  # [s]
const DTI = 1.0  # [s]
const G = 9.8  # [m/s^2]

const Xmin = -2000.0  # [m]
const Xmax = 2000.0  # [m]
const Ymin = -2000.0  # [m]
const Ymax = 2000.0  # [m]
const Pmin = 0.0  # [rad]
const Pmax = 2 * pi  # [rad]
const Vmin = 10  # [m/s]
const Vmax = 20  # [m/s]

const Xdim = 11
const Ydim = 11
const Pdim = 5
const Vdim = 3
const COCdim = 2

const NStates = Xdim * Ydim * Pdim * Vdim^2 * COCdim + 1
const NActions = 6^2

const Xs = linspace(Xmin, Xmax, Xdim)
const Ys = linspace(Ymin, Ymax, Ydim)
const Bearings = linspace(Pmin, Pmax, Pdim)
const Vs = linspace(Vmin, Vmax, Vdim)

end # module SCAConst