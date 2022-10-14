dATs = {2456109.5: 35, 2457204.5: 36, 2457754.5: 37}

def GetdAT(JD_UTC):
    # Retrieve all the dATs which start at JDs before the date
    # of interest. The last of these is the desired dAT.
    dAT_list = [v for k,v in dATs.items() if k < JD_UTC]
    
    return dAT_list[-1]
