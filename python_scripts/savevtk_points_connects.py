##############################################################################
#
# FUNCTION: prints unstructured grid (CONNECTED point) data to VTK formated file
#
#
# Author: Nicholas A. Battista
# Date: 8/24/16
# Github: http://github.org/nickabattista
# Institution: UNC-CH
# Lab: Laura Miller Lab
#
##############################################################################

def savevtk_points_connects( X, filename, vectorName,connectsMat):
    '''Prints matrix vector data to vtk formated file
    
    Args:
        X: Matrix of size Nx3
        filename: File name
        vectorname:
        connectsMat:'''

    # X is matrix of size Nx3 
    #              Col 1: x-data
    #              Col 2: y-data
    #              Col 3: z-data
    # filename:    What you are saving the VTK file as (string)
    # vectorName:  What you are naming the data you're printing (string)
    # connectsMat: Describes connections between nodes


    N = X.shape[0]
    Nc = connectsMat.shape[0]

    if C_flag==True:
        #Just add the measure of time for transforming the 
        nX = np.ascontiguousarray(X, dtype=np.float64)
        nconnectsMat = np.ascontiguousarray(connectsMat, dtype=np.float64)
        write.savevtk_points_connects_write(N,Nc,nX,filename,vectorName,nconnectsMat)

    else:
        with open(filename,'w') as file:
            file.write('# vtk DataFile Version 2.0\n')
            file.write(vectorName+'\n')
            file.write('ASCII\n')
            file.write('DATASET UNSTRUCTURED_GRID\n\n')
            #
            file.write('POINTS {0} float\n'.format(N))
            for ii in range(N):
                file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
            file.write('\n')
            #
            #First: # of "Cells", Second: Total # of info inputed following
            file.write('CELLS {0} {1}\n'.format(Nc,3*Nc))
            for s in range(Nc):
                file.write('{0} {1:d} {2:d}\n'.format(2,connectsMat[s,0],connectsMat[s,1]))
            file.write('\n')
            #
            file.write('CELL_TYPES {0}\n'.format(Nc)) # N = # of "Cells"
            for ii in range(Nc):
               file.write('3 ')
            file.write('\n')
    




