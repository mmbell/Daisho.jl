# Preliminary attempt at Daisho parameter list

# Raw moment names
raw_moment_names = ["DBZ", "ZDR", "KDP", "RHOHV", "VEL", "WIDTH", "PHIDP", "SQI", "SNR", "UNKNOWN_ID_82" ]

# QC moment names
qc_moment_names = ["DBZ", "ZDR", "KDP", "RHOHV", "VEL", "WIDTH", "PHIDP", "SQI"]
moment_grid_type = [:linear, :linear, :weighted, :weighted, :weighted, :weighted, :weighted, :weighted]
# Setting these as global variables now, but should be parameters
sqi_threshold = 0.35
snr_threshold = 6.0
spec_width_threshold = 8.0
rhohv_threshold = 0.7

# Regular grid parameters
# Eventually pass this to Springsteel.jl
xmin = -125000.0
xincr = 500.0
xdim = 501

ymin = -125000.0
yincr = 500.0
ydim = 501

zmin = 0.0
zincr = 500.0
zdim = 37

rmin = 0.0
rincr = 250.0
rdim = 481

rhi_zmin = 0.0
rhi_zincr = 250.0
rhi_zdim = 41

# Expand radius of influence along with the beam
beam_inflation = 0.01
power_threshold = 0.5
