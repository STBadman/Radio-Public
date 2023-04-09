import numpy as np
import astropy.units as u
import astropy.constants as const

def triangulate(dt,time_of_arrival_dict) :
    if not isinstance(time_of_arrival_dict,dict) :
        return "Input should be datetime instance and dictionary with three entries of spacecraft and times of arrival in seconds"
    sc_pos = list(locate(dt,time_of_arrival_dict.keys()).values())
    r01 = (np.array([sc_pos[1].x-sc_pos[0].x,
                    sc_pos[1].y-sc_pos[0].y]
                    )*sc_pos[0].x.unit).to("au").flatten()
    r02 = (np.array([sc_pos[2].x-sc_pos[0].x,
                    sc_pos[2].y-sc_pos[0].y]
                   )*sc_pos[0].x.unit).to("au").flatten()
    r12 = (np.array([sc_pos[2].x-sc_pos[1].x,
                    sc_pos[2].y-sc_pos[1].y]
                  )*sc_pos[0].x.unit).to("au").flatten()
    
    toas = list(time_of_arrival_dict.values())
    toas = [np.array(toa) if isinstance(toa,float) else toa for toa in toas]
    ct01 = ((toas[1]-toas[0])*const.c*u.s).to("au")
    ct02 = ((toas[2]-toas[0])*const.c*u.s).to("au")
    ct12 = ((toas[2]-toas[1])*const.c*u.s).to("au")

    # If measured delay of arrival larger than time to prop between sc return nan
    for jj,ct in enumerate(ct01) :
        if (ct/np.linalg.norm(r01)).value > 1 : ct01[jj]=np.nan
    for jj,ct in enumerate(ct02) :
        if (ct/np.linalg.norm(r02)).value > 1 : ct02[jj]=np.nan
    
    # Properties of baseline vectors in input frame.
    rmag01,rmag02,rmag12=np.linalg.norm(r01),np.linalg.norm(r02),np.linalg.norm(r12)
    rhat01,rhat02,rhat12 = r01/rmag01,r02/rmag02,r12/rmag12
    theta102 = np.arccos(np.dot(rhat01,rhat02))
    theta102 *= np.sign(rhat01[0]*rhat02[1] - rhat01[1]*rhat02[0])

    # Hyperbola parameterization values
    a_x = rmag01/2.*(1.-(ct01/rmag01).value**2)
    b_x = (ct01/rmag01).value
    a_y = (rmag02/2.*(1.-(ct02/rmag02).value**2)-a_x*np.cos(theta102))/np.sin(theta102)
    b_y = ((ct02/rmag02).value-b_x*np.cos(theta102))/np.sin(theta102)

    # Linear Parameterization of Source 
    r_WS1 = 1./2./(b_x**2+b_y**2-1)*(-2*(a_x*b_x+a_y*b_y) + (4.*(a_x*b_x+a_y*b_y
            )**2 - 4.*(b_x**2+b_y**2-1)*(a_x**2+a_y**2))**0.5)
    r_WS2 = 1./2./(b_x**2+b_y**2-1)*(-2*(a_x*b_x+a_y*b_y) - (4.*(a_x*b_x+a_y*b_y
            )**2 - 4.*(b_x**2+b_y**2-1)*(a_x**2+a_y**2))**0.5)  
    x1 = a_x + b_x*r_WS1
    y1 = a_y + b_y*r_WS1
    x2 = a_x + b_x*r_WS2
    y2 = a_y + b_y*r_WS2

    rshift = np.array([sc_pos[0].x.to("au"),sc_pos[0].y.to("au")])

    # Rotate Solution to Input Frame
    rot2hee = rot_2D(np.arctan2(rhat01[1],rhat01[0]))
    [x1,y1] = np.dot(rot2hee,[x1,y1])
    [x2,y2] = np.dot(rot2hee,[x2,y2])


    return  (
            (np.array([x1+rshift[0],y1+rshift[1]]))*u.au,
            (np.array([x2+rshift[0],y2+rshift[1]]))*u.au,
            ) 