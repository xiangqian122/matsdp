def input_default():
    dvm_input_default_dict = {}
    dvm_input_default_dict['jobid'] = '\"dvm\"'
    dvm_input_default_dict['unitscale']    = 1.889      
    dvm_input_default_dict['nbasis']       = 1          
    dvm_input_default_dict['natom']        = 1        
    dvm_input_default_dict['nsystem']      = 1          
    dvm_input_default_dict['basisfile']    = '\"D53.DAT\"'  
    dvm_input_default_dict['frozen']      = 1           
    dvm_input_default_dict['functional']   = 0          
    dvm_input_default_dict['spin']         = 0          
    dvm_input_default_dict['charge']       = 0.00       
    dvm_input_default_dict['ef']           = 0.00       
    dvm_input_default_dict['beta']         = 4000.000   
    dvm_input_default_dict['pftype']       = 0          
    dvm_input_default_dict['ipop']         = 1          
    dvm_input_default_dict['isolve']       = 1          
    dvm_input_default_dict['ifefld']       = 0          
    dvm_input_default_dict['efld']         = '0.0 0.0 0.0'
    dvm_input_default_dict['scfstep']      = 800 
    dvm_input_default_dict['mixtype']      = 0          
    dvm_input_default_dict['occshift']     = 0.0002     
    dvm_input_default_dict['ifscfdiis']    = 1          
    dvm_input_default_dict['nscfdiis']     = 12         
    dvm_input_default_dict['scfdiistol']   = 0.5        
    dvm_input_default_dict['scfrhotol']    = 0.00010    
    dvm_input_default_dict['scfspntol']    = 0.00010    
    dvm_input_default_dict['scfrhomix']    = 0.01       
    dvm_input_default_dict['scfspnmix']    = 0.01       
    dvm_input_default_dict['grid3dseed']   = '2 3 5'      
    dvm_input_default_dict['nbytediv']     = 1000       
    dvm_input_default_dict['ifbsfn']       = 0          
    dvm_input_default_dict['smplmixload']  = 0          
    dvm_input_default_dict['ifetot']       = 1          
    dvm_input_default_dict['ifbypscfc']    = 0          
    dvm_input_default_dict['printbasis']   = 0          
    dvm_input_default_dict['printscf']     = 3          
    dvm_input_default_dict['ifpopu']       = 1          
    dvm_input_default_dict['ifpdos']       = 0          
    dvm_input_default_dict['nplot']        = 0
    dvm_input_default_dict['plotinfo']     = 111        
    dvm_input_default_dict['ifdebug']      = 0          
    dvm_input_default_dict['frozeninfo']   = ''
    dvm_input_default_dict['grid3dinfo']   = '   600   1.00000   2.00000'
    dvm_input_default_dict['plotgrid']     ='''     2   100   100             &    ! (1) plot%dim, plot%ngrid
                264.600  47.5800  243.300       &
                264.600  47.5800  259.400       &
                247.300  47.5800  243.300       &
                   2   100   100                &    ! (2) plot%dim, plot%ngrid
                264.600  47.5800  243.300       &
                264.600  47.5800  227.200       &        
                247.300  47.5800  243.300       &
                   2   100   100                &    ! (3) plot%dim, plot%ngrid               
                264.600  47.5800  227.200       &
                264.600  47.5800  211.100       &                        
                247.300  47.5800  227.200       &
                   2   100   100                &    ! (4) plot%dim, plot%ngrid             
                264.600  47.5800  259.400       &
                264.600  47.5800  275.500       &                                  
                247.300  47.5800  259.400       &
                   2   100   100                &    ! (5) plot%dim, plot%ngrid
                264.600  47.5800  275.500       &
                264.600  47.5800  291.600       &
                247.300  47.5800  275.500'''
    return dvm_input_default_dict
