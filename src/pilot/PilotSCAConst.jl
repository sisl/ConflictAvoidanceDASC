module PilotSCAConst

export PenConflict, PenMinSep, PenCloseness, PenAction, StdDist, MinSepSq, InvVar
export DT, DTI, G, Xmin, Xmax, Ymin, Ymax, Bearingmin, Bearingmax, Speedmin, Speedmax
export Xdim, Ydim, Bearingdim, Speeddim, COCdim, NStates, NActions
export Xs, Ys, Bearings, Speeds, Actions
export SigmaSpeed, SigmaBank, SigmaBankCOC
export SigmaDim, SigmaWeightNominal, SigmaWeightOffNominal
export Responses, MeanResponseTime, RespondingProb, NonRespondingProb


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

const Xdim = 41  # 11
const Ydim = 41  # 11
const Bearingdim = 37  # 5
const Speeddim = 3

const Responsedim = 2

const NStates = Xdim * Ydim * Bearingdim * Speeddim^2 * Responsedim^2 + 1
const NActions = 36

const Xs = linspace(Xmin, Xmax, Xdim)
const Ys = linspace(Ymin, Ymax, Ydim)
const Bearings = linspace(Bearingmin, Bearingmax, Bearingdim)
const Speeds = linspace(Speedmin, Speedmax, Speeddim)

const Responses = [0.0, 1.0]  # indicator values: false == 0, true == 1

const MeanResponseTime = 5.0  # [s]
const RespondingProb = DT / (DT + MeanResponseTime)
const NonRespondingProb = 1.0 - RespondingProb

const Actions = [:right20, :right10, :straight, :left10, :left20, :clearOfConflict]

const SigmaSpeed = 2.0  # [m/s]
const SigmaBank = deg2rad(4.0)  # [rad]
const SigmaBankCOC = deg2rad(10.0)  # [rad]

const SigmaDim = 4  # number of dimensions for sigma-point sampling
const SigmaWeightNominal = 1 / 3
const SigmaWeightOffNominal = (1 - SigmaWeightNominal) / (2 * SigmaDim)

end # module PilotSCAConst