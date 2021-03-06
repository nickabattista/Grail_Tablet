##############################################################################
#
# FUNCTION: prints unstructured point data w/ associated scalar value to vtk 
#           formated file
#
#
# Author: Nicholas A. Battista
# Date: 8/24/16
# Github: http://github.org/nickabattista
# Institution: UNC-CH
# Lab: Laura Miller Lab
#
##############################################################################

def savevtk_points_with_scalar_data( X, scalarArray, filename, colorMap):
    ''' Prints matrix vector data to vtk formated file
    
    Args:
        X: Matrix of size Nx3
        scalarArray:
        filename:
        colorMap:'''

    # X is matrix of size Nx3 
    #              Col 1: x-data
    #              Col 2: y-data
    #              Col 3: z-data
    # scalarArray: Scalar array you are assigning to each point
    # filename:    What you are saving the VTK file as (string)
    # colorMap:  What you are naming the data you're printing (string)

    N = X.shape[0]
    nx = scalarArray.shape[0]

    if C_flag == True:
        nX = np.ascontiguousarray(X, dtype=np.float64)
        write.savevtk_points_write(N,nX,filename,colorMap)
    else:
        with open(filename,'w') as file:
            file.write('# vtk DataFile Version 2.0\n')
            file.write(colorMap+'\n')
            file.write('ASCII\n')
            file.write('DATASET UNSTRUCTURED_GRID\n\n')
            file.write('POINTS {0} float\n'.format(N))
            for ii in range(N):
                file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
            file.write('\n')
            #
            file.write('POINT_DATA   {0}\n'.format(nx*1*1))
            file.write('SCALARS '+colorMap+' double\n')
            file.write('LOOKUP_TABLE default\n')
            file.write('\n')
            for c in range(nx):
                file.write('{0} '.format(scalarArray[c]))
                file.write('\n')
            #Python 3.5 automatically opens in text mode unless otherwise specified