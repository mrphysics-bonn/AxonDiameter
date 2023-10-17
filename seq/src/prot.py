""" Helper functions for writing protocol data via Ismrmrd or (deprecated) h5py.
"""

import ismrmrd

def create_hdr(hdr, params):

    # initialize some optional parameters
    if "dwelltime" not in params:
        params["dwelltime"] = 0
    if "slices" not in params:
        params["slices"] = 1
    if "nsegments" not in params:
        params["nsegments"] = 1
    if "npartitions" not in params:
        params["npartitions"] = 1
    if "nintl" not in params:
        params["nintl"] = 1
    if "avg" not in params:
        params["avg"] = 1
    if "rep" not in params:
        params["rep"] = 1
    if "ncontrast" not in params:
        params["ncontrast"] = 1
    if "traj_delay" not in params:
        params["traj_delay"] = 0
    if "t_min" not in params:
        params["t_min"] = 0
    if "os_region" not in params:
        params["os_region"] = 1
    if "os_factor" not in params:
        params["os_factor"] = 1
    if "redfac" not in params:
        params["redfac"] = 1
    if "sms_factor" not in params:
        params["sms_factor"] = 1
    if "slice_res" not in params:
        params["slice_res"] = 1
    if "fov_z" not in params:
        params["fov_z"] = 1
    if "half_refscan" not in params:
        params["half_refscan"] = False
    if "trajtype" not in params:
        raise ValueError("Specify trajectory type: spiral or cartesian")
        
    # experimental conditions
    exp = ismrmrd.xsd.experimentalConditionsType()
    exp.H1resonanceFrequency_Hz = 298060000 # 7T
    hdr.experimentalConditions = exp

    # user parameters
    dtime = ismrmrd.xsd.userParameterDoubleType()
    dtime.name = 'dwellTime_us'
    dtime.value = params["dwelltime"] * 1e6 
    tdel = ismrmrd.xsd.userParameterDoubleType()
    tdel.name = 'traj_delay'
    tdel.value = params["traj_delay"] # this field can be used for any additional trajectory delay [s]
    nseg = ismrmrd.xsd.userParameterDoubleType()
    nseg.name = 'nsegments'
    nseg.value = params["nsegments"] # deprecated, but kept for compatibility
    tmin = ismrmrd.xsd.userParameterDoubleType()
    tmin.name = 't_min'
    tmin.value = params["t_min"] # starting time of t-vector for B0-correction
    os_region = ismrmrd.xsd.userParameterDoubleType()
    os_region.name = 'os_region'
    os_region.value = params["os_region"] # factor for the spiral oversampling region - phase maps will be reconstructed up to this factor of kspace
    os_factor = ismrmrd.xsd.userParameterDoubleType()
    os_factor.name = 'os_factor'
    os_factor.value = params["os_factor"] # ADC oversampling factor of spiral acquisition

    hdr.userParameters = ismrmrd.xsd.userParametersType()
    hdr.userParameters.userParameterDouble.append(dtime)
    hdr.userParameters.userParameterDouble.append(tdel)
    hdr.userParameters.userParameterDouble.append(nseg)
    hdr.userParameters.userParameterDouble.append(tmin)
    hdr.userParameters.userParameterDouble.append(os_region)
    hdr.userParameters.userParameterDouble.append(os_factor)

    # encoding
    encoding = ismrmrd.xsd.encodingType()
    if params["trajtype"] == "spiral":
        encoding.trajectory = ismrmrd.xsd.trajectoryType('spiral')
    if params["trajtype"] == "cartesian":
        encoding.trajectory = ismrmrd.xsd.trajectoryType('cartesian')

    # encoded and recon spaces
    efov = ismrmrd.xsd.fieldOfViewMm()
    efov.x = params["fov"]
    efov.y = params["fov"]
    efov.z = 1e3 * params["slice_res"] # used for image position calculation (should include the distance factor)
    rfov = ismrmrd.xsd.fieldOfViewMm()
    rfov.x = params["fov"]
    rfov.y = params["fov"]
    rfov.z = 1e3 * params["fov_z"] # used for scaling of SMS trajectories

    ematrix = ismrmrd.xsd.matrixSizeType()
    ematrix.x = int(params["fov"]/params["res"]+0.5)
    ematrix.y = int(params["fov"]/params["res"]+0.5)
    ematrix.z = params["npartitions"]
    rmatrix = ismrmrd.xsd.matrixSizeType()
    rmatrix.x = int(params["fov"]/params["res"]+0.5)
    rmatrix.y = int(params["fov"]/params["res"]+0.5)
    rmatrix.z = params["npartitions"]
    
    espace = ismrmrd.xsd.encodingSpaceType() # encoded space incl oversampling
    espace.matrixSize = ematrix
    espace.fieldOfView_mm = efov
    rspace = ismrmrd.xsd.encodingSpaceType()
    rspace.matrixSize = rmatrix
    rspace.fieldOfView_mm = rfov
    
    # set encoded and recon spaces
    encoding.encodedSpace = espace
    encoding.reconSpace = rspace
    
    # encoding limits
    limits = ismrmrd.xsd.encodingLimitsType()
    limits.slice = ismrmrd.xsd.limitType()
    limits.slice.minimum = 0
    limits.slice.maximum = params["slices"] - 1
    limits.slice.center = int(params["slices"]/2)

    limits.kspace_encoding_step_1 = ismrmrd.xsd.limitType()
    limits.kspace_encoding_step_1.minimum = 0
    limits.kspace_encoding_step_1.maximum = params["nintl"] - 1
    limits.kspace_encoding_step_1.center = 0

    limits.average = ismrmrd.xsd.limitType()
    limits.average.minimum = 0
    limits.average.maximum = params["avg"] - 1
    limits.average.center = 0

    limits.repetition = ismrmrd.xsd.limitType()
    limits.repetition.minimum = 0
    limits.repetition.maximum = params["rep"] - 1
    limits.repetition.center = 0

    limits.phase = ismrmrd.xsd.limitType()
    limits.phase.minimum = 0
    limits.phase.maximum =  0
    limits.phase.center = 0

    limits.contrast = ismrmrd.xsd.limitType()
    limits.contrast.minimum = 0
    limits.contrast.maximum = params["ncontrast"] - 1 # all contrasts that were acquired (relevant parameter for PowerGrid)
    limits.contrast.center = 0

    limits.segment = ismrmrd.xsd.limitType()
    limits.segment.minimum = 0
    limits.segment.maximum = params["nsegments"] - 1 # number of ADC segments
    if  params['half_refscan']:
        limits.segment.center = 1 # used as indicator for halved number of refscan slices
    else:
        limits.segment.center = 0

    pi = ismrmrd.xsd.parallelImagingType()
    acc_fac = ismrmrd.xsd.accelerationFactorType()
    acc_fac.kspace_encoding_step_1 = params['redfac']
    acc_fac.kspace_encoding_step_2 = params['sms_factor'] # use this field for multiband factor
    pi.accelerationFactor = acc_fac
    encoding.parallelImaging = pi

    encoding.encodingLimits = limits
    
    hdr.encoding.append(encoding)
