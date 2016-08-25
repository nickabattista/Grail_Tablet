##############################################################################
#
# FUNCTION: prints structured point vector data to vtk formated file
#
#
# Author: Nicholas A. Battista
# Date: 8/24/16
# Github: http://github.org/nickabattista
# Institution: UNC-CH
# Lab: Laura Miller Lab
#
##############################################################################

def savevtk_vector(X, Y, filename, vectorName,dx,dy):
    ''' Prints matrix vector data to vtk formated file.
    
    Args:
        X: 2-D ndarray
        Y: 2-D ndarray
        filename: file name
        vectorName:
        dx:
        dy:'''
    
    # X,Y is matrix of size Nx3 containing X,Y-direction vector components
    #              Col 1: x-data
    #              Col 2: y-data
    #              Col 3: z-data
    # filename:    What you are saving the VTK file as (string)
    # vectorName:  What you are naming the data you're printing (string)
    # dx,dy:       Grid spacing (resolution)

    #   Note:
    #   3D is clearly broken in this code, but there were still some reminants 
    #   in the matlab version. Given the choice of doing try/except blocks to
    #   keep these reminants or to kill them entirely, I'm choosing to kill them.
    #   So, specifically, nz is now gone. I will keep the output the same,
    #   however, for compatibility. So 1 will be printed in the Z column.
    
    #   Checks for compatibility of array sizes
    assert (X.shape == Y.shape), 'Error: velocity arrays of unequal size'
    nx, ny = X.shape
    
    XRow = X.shape[0]
    XCol = X.shape[1]
    YRow = Y.shape[0]
    YCol = Y.shape[1]


    if C_flag == True:
        nX = np.ascontiguousarray(X, dtype=np.float64)
        nY = np.ascontiguousarray(Y, dtype=np.float64)
        write.savevtk_vector(XRow,XCol,YRow,YCol,nX,nY,filename,vectorName,dx,dy)
    else:
        with open(filename,'w') as fid:
            fid.write('# vtk DataFile Version 2.0\n')
            fid.write('Comment goes here\n')
            fid.write('ASCII\n')
            fid.write('\n')
            fid.write('DATASET STRUCTURED_POINTS\n')
            # 1 below was nz
            fid.write('DIMENSIONS    {0}   {1}   {2}\n'.format(nx, ny, 1))
            fid.write('\n')
            fid.write('ORIGIN    0.000   0.000   0.000\n')
            #fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
            fid.write('SPACING   '+str(dx)+str(' ')+str(dy)+'   1.000\n')
            fid.write('\n')
            fid.write('POINT_DATA   {0}\n'.format(nx*ny))
            fid.write('VECTORS '+vectorName+' double\n')
            fid.write('\n')
            for b in range(ny):
                for c in range(nx):
                    fid.write('{0} '.format(X[c,b]))
                    fid.write('{0} '.format(Y[c,b]))
                    fid.write('0 ')
                fid.write('\n')
    #Python 3.5 automatically opens in text mode unless otherwise specified
