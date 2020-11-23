def wcp(V_step=-5, step_start=10, step_duration=20):
    """
    Measures whole cell properties. Specifically, this function returns the
    voltage clamp step estimates of series resistance, input resistance, cell 
    membrane resistance, cell membrane capacitance, cell surface area and 
    specific membrane resistance.
    
    The series (or access) resistance is obtained my dividing the voltage step
    by the peak amplitude of the current transient (Ogden, 1994): Rs = V / Ip
    
    The input resistance is obtained by dividing the voltage step by the average 
    amplitude of the steady-state current (Barbour, 2014): Rin = V / Iss
    
    The cell membrane resistance is calculated by subtracting the series 
    resistance from the input resistance (Barbour, 1994): Rm = Rin - Rs
    
    The cell membrane capacitance is estimated by dividing the transient charge 
    by the size of the voltage-clamp step (Taylor et al. 2012): Cm = Q / V
    
    The cell surface area is estimated by dividing the cell capacitance by the
    specific cell capacitance, c (1.0 uF/cm^2; Gentet et al. 2000; Niebur, 2008):
    Area = Cm / c
    
    The specific membrane resistance is calculated by multiplying the cell
    membrane resistance with the cell surface area: rho = Rm * Area
    Users should be aware of the approximate nature of determining cell
    capacitance and derived parameters from the voltage-clamp step method
    (Golowasch, J. et al., 2009)
    References:
    Barbour, B. (2014) Electronics for electrophysiologists. Microelectrode
     Techniques workshop tutorial.
     www.biologie.ens.fr/~barbour/electronics_for_electrophysiologists.pdf
    Gentet, L.J., Stuart, G.J., and Clements, J.D. (2000) Direct measurement
     of specific membrane capacitance in neurons. Biophys J. 79(1):314-320
    Golowasch, J. et al. (2009) Membrane Capacitance Measurements Revisited: 
     Dependence of Capacitance Value on Measurement Method in Nonisopotential 
     Neurons. J Neurophysiol. 2009 Oct; 102(4): 2161-2175.
    Niebur, E. (2008), Scholarpedia, 3(6):7166. doi:10.4249/scholarpedia.7166
     www.scholarpedia.org/article/Electrical_properties_of_cell_membranes
     (revision #13938, last accessed 30 April 2018)
    Ogden, D. Chapter 16: Microelectrode electronics, in Ogden, D. (ed.) 
     Microelectrode Techniques. 1994. 2nd Edition. Cambridge: The Company
     of Biologists Limited.
    Taylor, A.L. (2012) What we talk about when we talk about capacitance 
     measured with the voltage-clamp step method J Comput Neurosci. 
     32(1):167-175
    """

    # Error checking
    if stf.get_yunits() != "pA":
        raise ValueError('The recording is not voltage clamp')
    
    # Prepare variables from input arguments
    si = stf.get_sampling_interval()
    t0 = step_start / si
    l = step_duration / si
    
    # Set cursors and update measurements
    stf.set_base_start((step_start - 1) / si)
    stf.set_base_end(t0-1)
    stf.set_peak_start(t0)
    stf.set_peak_end((step_start + 1) / si)
    stf.set_fit_start(t0)
    stf.set_fit_end(t0+l-1)
    stf.set_peak_direction("both")
    stf.measure()

    # Calculate series resistance (Rs) from initial transient
    b = stf.get_base()
    Rs = 1000 * V_step / (stf.get_peak() - b)  # in Mohm

    # Calculate charge delivered during the voltage clamp step
    n = int(stf.get_fit_end()+1-stf.get_fit_start())
    x = [i*stf.get_sampling_interval() for i in range(n)]
    y = stf.get_trace()[int(stf.get_fit_start()):int(stf.get_fit_end()+1)]
    Q = np.trapz(y-b,x)

    # Set cursors and update measurements
    stf.set_base_start(t0+l-1-(step_duration/4)/si)
    stf.set_base_end(t0+l-1)
    stf.measure()

    # Measure steady state current and calculate input resistance
    I = stf.get_base() - b
    Rin  = 1000 * V_step / I                   # in Mohm
    
    # Calculate cell membrane resistance   
    Rm = Rin - Rs                              # in Mohm

    # Calculate voltage-clamp step estimate of the cell capacitance
    t = x[-1] - x[0]
    Cm = (Q - I * t) / V_step                  # in pF

    # Estimate membrane surface area, where the capacitance per unit area is 1.0 uF/cm^2
    A = Cm * 1e-06 / 1.0                       # in cm^2
    
    # Calculate specific membrane resistance
    rho = 1e+03 * Rm * A                       # in kohm.cm^2; usually 10 at rest

    # Create table of results
    retval  = []
    retval += [("Holding current (pA)", b)]
    retval += [("Series resistance (Mohm)", Rs)]
    retval += [("Input resistance (Mohm)", Rin)]
    retval += [("Cell resistance (Mohm)", Rm)]
    retval += [("Cell capacitance (pF)", Cm)]
    retval += [("Surface area (um^2)", A * 1e+04**2)]
    retval += [("Membrane resistivity (kohm.cm^2)", rho)]
    retval = dict(retval)
    stf.show_table(retval,"Whole-cell properties")
    
    return
