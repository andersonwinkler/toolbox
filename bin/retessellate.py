#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:11:57 2024

@author: winkleram
"""

import sys
import os
import nibabel as nib
import numpy as np

if __name__ == "__main__":

    # Ensure correct number of arguments
    if len(sys.argv) != 5:
        print("Usage: script.py <path1> <path2> <path3> <path4>")
        sys.exit(1)

    # Unpack arguments
    path1, path2, path3, path4 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

    # Default margin
    marg =  0.05;
    
    # Check if the first three paths exist.
    if not os.path.exists(path1):
        raise FileNotFoundError(f"Path does not exist: {path1}")
    if not os.path.exists(path2):
        raise FileNotFoundError(f"Path does not exist: {path2}")
    if not os.path.exists(path3):
        raise FileNotFoundError(f"Path does not exist: {path3}")

    # Load the surfaces using nibabel for FreeSurfer formats
    vtx1, fac1 = nib.freesurfer.read_geometry(path1)
    vtx2, fac2 = nib.freesurfer.read_geometry(path2)
    vtx3, fac3 = nib.freesurfer.read_geometry(path3)
    nF1 = fac1.shape[0]
    nV2 = vtx2.shape[0]
    
    # Sanity check
    if not np.array_equal(fac1, fac3):
        raise ValueError("The source surface and the surface to be retessellated must have the same geometry.")

    # Where the result is going to be stored
    vtx4 = np.zeros((nV2, 3))
    
    # Vertices' coords per face
    facvtx1 = np.hstack((vtx1[fac1[:, 0], :], vtx1[fac1[:, 1], :], vtx1[fac1[:, 2], :]))
    
    # Face barycenter
    xbary = np.mean(facvtx1[:, [0, 3, 6], None], axis=1)    # x-coordinate
    ybary = np.mean(facvtx1[:, [1, 4, 7], None], axis=1)    # y-coordinate
    zbary = np.mean(facvtx1[:, [2, 5, 8], None], axis=1)    # z-coordinate
    cbary = np.hstack((xbary, ybary, zbary))                # Cartesian coordinates of the barycenters
    r     = np.sqrt(xbary**2 + ybary**2 + zbary**2)         # radius
    theta = np.arctan2(ybary, xbary)                        # azimuth (angle in x-y plane)
    phi   = np.arctan2(zbary, np.sqrt(xbary**2 + ybary**2)) # elevation (angle from z-axis)
    sbary = np.hstack((theta, phi, r))                      # Spherical coordinates of the barycenters
    
    # Pre-calculated sines and cosines of azimuth and elevation:
    sinA = np.sin(sbary[:, 0, None])
    sinE = np.sin(sbary[:, 1, None])
    cosA = np.cos(sbary[:, 0, None])
    cosE = np.cos(sbary[:, 1, None])
    
    # Pre-calculated rotation matrices
    rotM = np.column_stack((cosA * cosE, sinA*cosE, sinE, -sinA, cosA, np.zeros(nF1), -cosA*sinE, -sinA*sinE, cosE))
    
    # Random angle around X
    rndangX = np.random.rand() * np.pi
    rndangX = np.pi/3
    sinX = np.sin(rndangX)
    cosX = np.cos(rndangX)
    rotM = np.column_stack((
        rotM[:, :3],  
        rotM[:, 3] * cosX + rotM[:, 6] * sinX,
        rotM[:, 4] * cosX + rotM[:, 7] * sinX,
        rotM[:, 5] * cosX + rotM[:, 8] * sinX,
        rotM[:, 6] * cosX - rotM[:, 3] * sinX,
        rotM[:, 7] * cosX - rotM[:, 4] * sinX,
        rotM[:, 8] * cosX - rotM[:, 5] * sinX 
    ))
    
    # Pre-calculated min and max for each face and bounding box
    minF = np.column_stack(( np.min(facvtx1[:, [0, 3, 6]], axis=1), np.min(facvtx1[:, [1, 4, 7]], axis=1), np.min(facvtx1[:, [2, 5, 8]], axis=1) ))
    maxF = np.column_stack(( np.max(facvtx1[:, [0, 3, 6]], axis=1), np.max(facvtx1[:, [1, 4, 7]], axis=1), np.max(facvtx1[:, [2, 5, 8]], axis=1) ))
    b    = np.tile(np.max((maxF - minF), axis=1).reshape(-1, 1), (1, 3)) * marg
    minF = minF - b
    maxF = maxF + b

    # For each source face
    for f in range(nF1):
    
        vidx = fac1[f, :]         # Indices of the vertices for face f
        Fvtx = vtx1[vidx, :]      # Corresponding vertex coordinates from vtx1
        
        # Candidate vertices
        Cidx = np.all((vtx2 >= minF[f, :]) & (vtx2 <= maxF[f, :]), axis=1)  # Logical condition across columns
        Cvtx = vtx2[Cidx, :]  # Extract candidate vertices
        Cidxi = np.where(Cidx)[0]  # Indices of the candidate vertices
        
        # Concatenate the face vertices and candidate vertices
        Avtx = np.vstack((Fvtx, Cvtx)) @ rotM[f, :].reshape(3, 3).T
        
        # Here is the main difference in relation to the 'distributive' method.
        # Instead of extrapolate (split) a point in the source into the face
        # vertices of the target for increments, find in the source
        # face what are the Target vertices that lie inside it and do a
        # barycentric interpolation (not to be confused with the rotation
        # of the barycenter, used some lines above to put the face under
        # analysis near the sphere equator and the meridian zero).
        
        # Convert to azimuthal gnomonic
        Gvtx = np.ones_like(Avtx)  # The 3rd col will remain full of ones
        Gvtx[:, 0] = Avtx[:, 1] / Avtx[:, 0]  # Tangent of the angle on the XY plane
        Gvtx[:, 1] = Avtx[:, 2] / Avtx[:, 0]  # Tangent of the angle on the XZ plane
        T = Gvtx[:3, :] # Face coords for the test below
        aT = np.linalg.det(T) # Face total area (2x the area, actually)
        
        # For every candidate vertex
        for v in range(len(Cidxi)):
            
            # Compute the areas for the subtriangles (2x area actually)
            # Subtriangle A
            tA = T.copy()  # Copy the original T
            tA[0, :] = Gvtx[v + 3, :]
            aA = abs(np.linalg.det(tA))
            
            # Subtriangle B
            tB = T.copy()
            tB[1, :] = Gvtx[v + 3, :]
            aB = abs(np.linalg.det(tB))
            
            # Subtriangle C
            tC = T.copy()
            tC[2, :] = Gvtx[v + 3, :]
            aC = abs(np.linalg.det(tC))
            
            # Test if the point is inside the face
            if np.float32(aT) == np.float32(aA + aB + aC): # Use float32 (single) to emulate Matlab. However, np.isclose would have been more pythonic
            
                # Weight by the areas and interpolate the value between the 3 vertices
                vtx4[Cidxi[v], :] = np.dot([aA, aB, aC], vtx3[vidx, :]) / aT
        
    # Save output surface
    nib.freesurfer.write_geometry(path4, vtx4, fac2)
