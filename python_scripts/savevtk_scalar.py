##############################################################################
#
# FUNCTION: prints 3D scalar matrix to VTK formated file
#
#
# Author: Nicholas A. Battista
# Date: 8/24/16
# Github: http://github.org/nickabattista
# Institution: UNC-CH
# Lab: Laura Miller Lab
#
##############################################################################

def savevtk_scalar(array, filename, colorMap,dx,dy):
    ''' Prints scalar matrix to vtk formatted file.
    
    Args:
        array: 2-D ndarray
        filename: file name
        colorMap:
        dx:
        dy:'''
    
    # array is matrix of size 'nx x ny x nz' containing scalar data on
    #              computational grid
    # filename:    What you are saving the VTK file as (string)
    # colorMap:    What you are naming the data you're printing (string)
    # dx,dy:       Grid spacing (resolution)

    #  	Note:
    #   3D is clearly broken in this code, but there were still some reminants 
    #   in the matlab version. Given the choice of doing try/except blocks to
    #   keep these reminants or to kill them entirely, I'm choosing to kill them.
    #   So, specifically, nz is now gone. I will keep the output the same,
    #   however, for compatibility. So 1 will be printed in the Z column.

    nx,ny = array.shape
    if C_flag == True:
        narray = np.ascontiguousarray(array, dtype=np.float64)
        write.savevtk_scalar(nx,ny,narray,filename,colorMap,dx,dy)
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
            fid.write('SPACING   '+str(dx)+str(' ')+str(dy)+'   1.000\n')
            fid.write('\n')
            # The 1 below was nz
            fid.write('POINT_DATA   {0}\n'.format(nx*ny*1))
            fid.write('SCALARS '+colorMap+' double\n')
            fid.write('LOOKUP_TABLE default\n')
            fid.write('\n')
            for b in range(ny):
                for c in range(nx):
                    fid.write('{0} '.format(array[c,b]))
                fid.write('\n')
        #Python 3.5 automatically opens in text mode unless otherwise specified
