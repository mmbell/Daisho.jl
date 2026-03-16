# Preliminary attempt at Daisho parameter list

# Raw moment names
const raw_moment_names = ["DBZ", "ZDR", "KDP", "RHOHV", "VEL", "WIDTH", "PHIDP", "SQI", "SNR", "UNKNOWN_ID_82" ]

# QC moment names
const qc_moment_names = ["DBZ", "ZDR", "KDP", "RHOHV", "VEL", "WIDTH", "PHIDP", "SQI"]
const moment_grid_type = [:linear, :linear, :weighted, :weighted, :weighted, :weighted, :weighted, :weighted]
# Setting these as global constants now
const sqi_threshold = 0.35
const snr_threshold = 6.0
const spec_width_threshold = 8.0
const rhohv_threshold = 0.7

# Regular grid parameters
# Eventually pass this to Springsteel.jl
const xmin = -125000.0
const xincr = 500.0
const xdim = 501

const ymin = -125000.0
const yincr = 500.0
const ydim = 501

const zmin = 0.0
const zincr = 500.0
const zdim = 37

const rmin = 0.0
const rincr = 250.0
const rdim = 481

const rhi_zmin = 0.0
const rhi_zincr = 250.0
const rhi_zdim = 41

# Expand radius of influence along with the beam
const beam_inflation = 0.01
const power_threshold = 0.5
