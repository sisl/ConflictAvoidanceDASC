module SCAConst

export PenConflict, PenMinSep, PenCloseness, PenAction, StdDist, MinSepSq, InvVar
export DT, DTI, G, Xmin, Xmax, Ymin, Ymax, Bearingmin, Bearingmax, Speedmin, Speedmax
export Xdim, Ydim, Bearingdim, Speeddim, COCdim, NStates, NActions
export Xs, Ys, Bearings, Speeds, Actions

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
const Bearingmin = 0.0  # [rad]
const Bearingmax = 2 * pi  # [rad]
const Speedmin = 10  # [m/s]
const Speedmax = 20  # [m/s]

const Xdim = 11
const Ydim = 11
const Bearingdim = 5
const Speeddim = 3

const NStates = Xdim * Ydim * Bearingdim * Speeddim^2 + 1
const NActions = 6^2

const Xs = linspace(Xmin, Xmax, Xdim)
const Ys = linspace(Ymin, Ymax, Ydim)
const Bearings = linspace(Bearingmin, Bearingmax, Bearingdim)
const Speeds = linspace(Speedmin, Speedmax, Speeddim)

const Actions = [:right20, :right10, :straight, :left10, :left20, :clearOfConflict]

end # module SCAConst