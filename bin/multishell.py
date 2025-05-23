#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 00:28:38 2025

@author: winkleram
"""

import os
import sys
import re
import numpy as np
from scipy import optimize as scopt

# ######[ HELP TEXT ]##########################################################

def print_help():
    help_text = '''
Multishell Diffusion MRI Sampling Tool

This command generates, manipulates, converts formats, and displays multishell
diffusion MRI sampling schemes. It can create new optimal schemes, load existing
ones from various formats, add or remove b=0 volumes, reorder directions, and
save the schemes in different formats. It also plots the directions.

A new scheme can be created using a simultaneous optimization (--new simultaneous)
of all directions across all shells (Equations 2 and 4 of the paper), or by
incrementally optimizing (--new incremental) by adding one direction at a time
(page 1537 of the paper, 1st column, 2nd paragraph, just before the Results section).

The simultaneous optimization is superior, as can be demonstrated by plotting
the directions (--plot interactive). But the incremental optimization has the
advantage that, if the acquisition during the MRI session is interrupted, the
result may still be usable as the incremental scheme guarantees a reasonably
uniform coverage with any number of direction.

But why tradeoff? Here, we can do simultaneous optimization, then reorder 
(with the option --reorder) the directions optimally, to obtain roughly the
same benefit of the incremental scheme, which is a reasonably unform coverage
if the MRI session is interrupted.

If we already have a diffusion sampling scheme in FSL/BIDS format (bvec/bval),
MRtrix3 format (*.b), Siemens format (*.dvs), GE format (*.dat), Philips format
(*.txt), or in the format produced by Ed Caruyer's webtool (*.txt), we can load
these, manipulate them by adding b0 at the beginning, end, or interspersed
among the other directions (--addb0), remove b0 (--removeb0), and specify or
change the bvalues (--bvalues).

For plotting (--plot), these can opened interactively, rotated, and saved
manually from the plotting window, or simply saved to a file in a specified
format such as *.pdf or *.png without any interaction (good for scripting).

Usage: multishell [options]

Options:

--help
    Display this help message.

--new <method> <Ks> | <K> <S> <distr>
    Create a new set of directions using the specified method ('simultaneous'
    or 'incremental'). Provide either:
    - <Ks>: a list of integers, e.g., "[10,20,30]", specifying directions per shell.
    - <K>: total directions, <S>: number of shells, <distr>: 'uniform', 'linear', or 'quadratic'.
    The <distr> parameter specifies if the total number of shells should be
    distributed uniformly across shells, or linearly as a function of the shell
    number, or quadratically as a function of the shell number.
    Examples:
    --new simultaneous '[10,20,30]'
    --new incremental 60 3 linear

--alpha <value>
    Set alpha (0 to 1) for balancing intra/inter-shell repulsion forces when
    creating new multishell schemes. See Equation 5 in the paper.
    Default: 0.75 (this is different from the paper, which suggests 0.50)

--load <filename> <format>
    Load scheme from file in specified format ('fsl', 'mrtrix3', 'siemens',
    'ge', 'philips', 'caruyer'). For 'siemens' and 'ge', provide --bmax.

--reorder
    Optimally reorder directions.

--bvalues <list>
    Specify b-values for shells, e.g., '[1000,2000,3000]'. This is needed if
    creating new schemes and saving to a file format that will be used in the
    scanner. Do not include here bvalues=0; use --addb0 for that.

--addb0 <list>
    Add b=0 volumes: [start, interspersed, end], e.g., "[4,5,3]".
    For the interspersed (middle) block, the first b-value will always be 0,
    so if you also indicate a number to be added at the start, you will find
    that there is one extra there; to discount that, subtract one from the start.

--removeb0
    Remove existing b=0 volumes from a file that includes that information.

--bmax <value>
    Specify maximum b-value for formats that scale vectors (Siemens, GE) as
    opposed to storing b-values explicitly (Philips, FSL, MRtrix3)

--save <prefix> [format]
    Save scheme with given prefix in specified format (default: "caruyer").
    For a 'philips', for the file to be valid, the first direction must have
    b-value = 0, and further, no directions can be repeated, even for multiple
    b-values = 0 (so, the direction 0,0,0 cannot be used). The command will
    take care of that, using random directions even for b-values = 0.

--plot [format]
    Plot directions. If no format is specified, opens interactively a new window;
    otherwise, saves in "pdf", "png", "webp", etc.

Examples:
    
multishell --new incremental 60 3 uniform --bvalues "[1000,2000,3000]" \\
           --addb0 "[4,5,5]" --save my_scheme siemens
              
multishell --load my_input.dvs siemens --bmax 3000 \\
           --removeb0 --save my_output ge
              
multishell --load my_input fsl --save my_output philips

multishell --new simultaneous ["10,20,30]" --reorder \\
           --bvalues "[1000,2000,3000]" --plot

Supported formats:
- fsl: FSL bvec/bval files, also BIDS as produced by Chris Rorden's dcm2niix.
- mrtrix3: MRtrix3 gradient table.
- siemens: Siemens DVS (requires --bmax).
- ge: GE tensor.dat (requires --bmax).
- philips: Philips text format (unique directions).
- caruyer: Caruyer's webtool format.

Theis commands implements the method originally proposed by:
    
* Caruyer E, Lenglet C, Sapiro G, Deriche R. Design of multishell sampling
  schemes with uniform coverage in diffusion MRI. Magn Reson Med.
  2013 Jun;69(6):1534-40. doi: 10.1002/mrm.24736.
  Epub 2013 Apr 26. PMID: 23625329; PMCID: PMC5381389.

If you use this command in your research, please cite the original paper.

_____________________________________
Anderson M. Winkler
Univ. of Texas Rio Grande Valley
April/2025
http://brainder.org
    '''
    print(help_text)


# ######[ FILE PARSING ]#######################################################

def read_fsl(filename): # =====================================================
    '''
    Read FSL bvec/bval files, also produced by dcm2niix and adopted by BIDS.
    '''
    basename = os.path.splitext(filename)[0]
    bvecpath = f'{basename}.bvec'
    bvalpath = f'{basename}.bval'
    if os.path.exists(bvecpath):
        bvec = np.loadtxt(bvecpath).T
    if os.path.exists(bvalpath):
        bval = np.loadtxt(bvalpath)[:,None]
    return bvec, bval

def write_fsl(filename, bvec, bval): # ----------------------------------------
    '''
    Write FSL bvec/bval files.
    '''
    basename = os.path.splitext(filename)[0]
    bvecpath = f'{basename}.bvec'
    bvalpath = f'{basename}.bval'
    np.savetxt(bvecpath, bvec.T, fmt='%.6f')
    np.savetxt(bvalpath, bval.T, fmt='%d')

def read_mrtrix3(filename): # =================================================
    '''
    Read an MRtrix3 gradient table. Note this is similar to Philips format,
    with less restrictions (no uniqueness constraint).
    '''
    gradv = []
    gradb = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip().split()
            v = list(map(float, line[:3]))
            b = float(line[3])
            gradv.append(v)
            gradb.append(b)
    gradv = np.array(gradv)
    gradb = np.array(gradb)[:,None]
    return gradv, gradb

def write_mrtrix3(filename, gradv, gradb): # ----------------------------------
    '''
    Write an MRtrix3 gradient table.
    '''   
    with open(filename, 'w') as f:
        for v, b in zip(gradv, gradb):
            line = f'{v[0]:.7g} {v[1]:.7g} {v[2]:.7g} {b.item():g}\n'
            f.write(line)

def read_siemens(filename): # =================================================
    '''
    Read direction scheme from a Siemens DVS file.
    '''
    metadata = {}
    vectors  = []
    vpattern = re.compile(r'Vector\[\d+\]\s*=\s*\(([-\d\s.,]+)\)')
    mpattern = re.compile(r'^\s*([a-zA-Z]+)\s*=\s*([^\s#]+)')
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            mmatch = mpattern.match(line)
            if mmatch:
                key   = mmatch.group(1).strip()
                value = mmatch.group(2).strip()
                metadata[key] = value
                continue
            vmatch = vpattern.search(line)
            if vmatch:
                vals = vmatch.group(1).replace(',', ' ').split()
                vec = [float(n) for n in vals]
                vectors.append(vec)
    vectors = np.array(vectors)
    return vectors, metadata

def write_siemens(filename, vectors, metadata=None): # -----------------------------
    '''
    Write a direction scheme in Siemens DVS format.
    '''
    if metadata is None:
        metadata = {'CoordinateSystem': 'xyz', 'Normalization': 'none'}
    with open(filename, 'w') as f:
        f.write('[directions={}]\n'.format(len(vectors)))
        for key, value in metadata.items():
            f.write(f'{key} = {value}\n')
        for i, v in enumerate(vectors):
            f.write(f'Vector[{i}] = ({v[0]}, {v[1]}, {v[2]})\n')

def read_ge(filename): # ======================================================
    '''
    Read a GE direction scheme from tensor.dat format.
    '''
    with open(filename) as f:
        lines = [line.strip() for line in f if not line.strip().isdigit()]
    vectors = np.array([list(map(float, line.split())) for line in lines])
    return vectors

def write_ge(filename, vectors): # --------------------------------------------
    '''
    Write a GE direction scheme in tensor.dat format.
    '''
    with open(filename, 'w') as f:
        f.write(f"{len(vectors)}\n")
        for v in vectors:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")

def read_philips(filename): # =================================================
    '''
    Read a Philips direction scheme from text format.
    '''
    vectors = []
    bvalues = []
    with open(filename, 'r') as f:
        comment = f.readline()
        line = comment.strip().split()
        if  len(line) == 4 and \
            all(part.replace('.', '', 1).replace('-', '', 1).isdigit() for part in line):
            vector = list(map(float, line[:3]))
            bvalue = int(line[3]) # keep separate if we need to change to int in the future
            vectors.append(vector)
            bvalues.append(bvalue)
            comment = None
        for line in f:
            line   = line.strip().split()
            vector = list(map(float, line[:3]))
            bvalue = float(line[3]) # keep separate if we need to change to int in the future
            vectors.append(vector)
            bvalues.append(bvalue)
    vectors = np.array(vectors)
    bvalues = np.array(bvalues)[:,None]
    return vectors, bvalues, comment

def write_philips(filename, vectors, bvalues, comment=None): # ----------------
    '''
    Write a Philips direction scheme to text format.
    Duplicate vectors are never allowed. For bvalues=0, use any random vector as
    long as it's unique (not necessarily (0,0,0).
    '''
    if bvalues[0] != 0:
        raise ValueError('For Philips, the first b-value must be 0.')
    if np.unique(vectors, axis=0).shape[0] != vectors.shape[0]:
        raise ValueError('For Philips, duplicate vectors are not allowed, not even for b=0.')
    with open(filename, 'w') as f:
        if comment:
            f.write(comment + '\n')
        for vector, bvalue in zip(vectors, bvalues):
            f.write(f'{vector[0]:.3f}\t{vector[1]:.3f}\t{vector[2]:.3f}\t{int(bvalue.item())}\n')

def read_caruyer(filename): # =================================================
    '''
    Read a text file specifying optimal directions as generated by the webtool
    in Emmanuel Caruyer's website:
    http://www.emmanuelcaruyer.com/q-space-sampling.php
    '''
    shells  = []
    vectors = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            shells.append(int(parts[0]))  # shell number
            vectors.append([float(parts[1]), float(parts[2]), float(parts[3])])
    vectors = np.array(vectors)
    shells  = np.array(shells)[:,None]
    shells = shells - shells.min() # we want the shells in python to start from 0
    return vectors, shells

def write_caruyer(filename, vectors, shells): # -------------------------------
    '''
    Write directions in the format used by the webtool in EC's website,
    without the disclaimer at the top.
    '''
    shells = shells - shells.min() + 1 # indices in the file start from 1
    with open(filename, 'w') as f:
        f.write('#shell\tu_x\tu_y\tu_z\n')
        for shell, vector in zip(shells, vectors):
            f.write(f'{shell[0]}\t{vector[0]}\t{vector[1]}\t{vector[2]}\n')


# ######[ SIMULTANEOUS OPTIMIZATION ]##########################################

def simultaneous_optimization(Ks, alpha=0.5, maxiter=500, epsilon=1e-9): # ===
    '''
    Simultaneously optimize sampling directions on multiple shells using
    the method described in:
    * Caruyer E, Lenglet C, Sapiro G, Deriche R. Design of multishell sampling
      schemes with uniform coverage in diffusion MRI.
      Magnetic Resonance in Med. 2013 Jun;69(6):1534–1540.
    '''
    Ks        = np.array(Ks)
    K         = np.sum(Ks)
    W, shells = calc_weights(Ks, alpha) # precompute weights
    vectors   = random_vectors(K)
    svectors  = vectors.ravel(order='C') # they need be stacked for optimization
    result    = scopt.minimize(
        cost_all_vectors,
        svectors,
        args        = (K,W,),
        method      = 'SLSQP',
        jac         = grad_all_vectors,
        constraints = [{'type'   : 'eq',
                        'fun'    : eq_constraints,
                        'args'   : (K,)}],
        options     =  {'maxiter': maxiter,
                        'ftol'   : epsilon,
                        'disp'   : False})
    svectors = result.x
    vectors  = svectors.reshape((K,3))
    vectors  = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors, shells

def calc_weights(Ks, alpha): # ----------------------------------------
    """
    Compute a K by K symmetric weight matrix that encapsulates the various
    constants from Equations 2 and 4 from the paper.
    """
    # This function could be simpler, without iterating twice over range(S)
    # but it's kept as this for clarity and symmetry with the paper.
    K      = np.sum(Ks) # this can be inferred from 
    S      = len(Ks)
    shells = np.concatenate([np.full(Ks[s],s) for s in range(S)]) # starts at 0
    W = np.zeros((K, K))
    
    # V1: repulsion within shell, Eqn 2
    for s in range(S):
        idxs   = (shells == s)
        intraW = alpha / (S * Ks[s]**2)
        idxss  = np.outer(idxs,idxs) - np.diag(idxs.astype(float))
        W     += idxss * intraW
    
    # V2: repulsion vetween shells, Eqn 4
    interW = (1 - alpha) / K**2
    for s in range(S):
        idxs = (shells == s)
        for t in range(S):
            if s != t:
                idxt  = (shells == t)
                idxst = np.outer(idxs,idxt)
                idxts = np.outer(idxt,idxs)
                W += (idxst + idxts) * interW
    shells = shells[:,None] # shells start counting from 0
    return W, shells

def cost_all_vectors(svectors, K, W, epsilon=1e-9): # -------------------------
    '''
    Cost function (simultaneous optimization of all directions).
    Note that svectors are the stacked coordinates of the vectors.
    '''
    vectors = svectors.reshape((K,3))
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    # Equation 3 of the paper. For unit vectors, note that:
    #    ||u - v||^2  =  2 - 2.u.v
    #    ||u + v||^2  =  2 + 2.u.v
    D = np.dot(vectors, vectors.T)   
    energy = (1 / (2 + 2*D + epsilon)) + (1 / (2 - 2*D + epsilon))
    Energy = np.sum(W * energy)
    return Energy

def grad_all_vectors(svectors, K, W, epsilon=1e-9): # -------------------------
    '''
    Gradient of the cost function.
    '''
    vectors = svectors.reshape((K,3))
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    num1    = vectors[:,None,:] - np.transpose(vectors[:,None,:], axes=(1,0,2))
    num2    = vectors[:,None,:] + np.transpose(vectors[:,None,:], axes=(1,0,2))
    D       = np.dot(vectors, vectors.T)
    den1    = (2 - 2*D)**2 + epsilon # ||u - u||^4
    den2    = (2 + 2*D)**2 + epsilon # ||u + u||^4
    pderiv  = num1/den1[:,:,None] + num2/den2[:,:,None] # Eqn. 7 of the paper, partial derivatives for the 3 axes
    grad    = -2*np.sum(W[:,:,None]*pderiv, axis=1) # sum over axis 0 or 1 but not both (either u-v or v-u, but not both)
    grad    = grad.ravel(order='C')
    return grad

def eq_constraints(svectors, K, *args): # -------------------------------------
    '''
    Equality constraint (vector norm equals 1).
    '''
    vectors = svectors.reshape((K,3), order='C')
    eq      = np.linalg.norm(vectors, axis=1) - 1.0
    return  eq

def optimal_reordering(vectors, shells): # ====================================
    '''
    After we've done the simultaneous optimization, we may want to reorder the
    vectors and shells such that they are in an optimal sequence as if they
    had been generated incrementally. This is what this function does,
    therefore combining the best of both possibilities.
    '''
    _, Ks   = np.unique(shells, return_counts=True)
    S       = len(Ks)
    sphvec  = cart2sph(vectors)[:,1:] # theta and phi only
    idx     = [] # assigned directions
    seq     = list(range(len(vectors))) # unassigned directions
    while len(seq):
        Ksc    = np.array([np.sum(shells[idx] == s) for s in range(S)])
        shell  = next_shell(Ksc, Ks)
        energy = np.full((len(seq),1), np.inf)
        for i, sq in enumerate(seq):
            if shells[sq] == shell:
                energy[i] = cost_one_vector(sphvec[sq], shells[sq], sphvec[idx], shells[idx], np.ones((S,S)))
        sidx   = np.argmin(energy)
        idx.append(seq[sidx])
        seq.remove(seq[sidx])
    vectors = vectors[idx]
    shells  = shells[idx]
    return vectors, shells, idx


# ######[ INCREMENTAL OPTIMIZATION ]###########################################

def incremental_optimization(Ks, alpha=0.5, maxiter=500): # ===================
    '''
    Incrementally optimize sampling directions on multiple shells using the
    method described briefly in page 1537 (last paragraph before Results) of:
    * Caruyer E, Lenglet C, Sapiro G, Deriche R. Design of multishell sampling
      schemes with uniform coverage in diffusion MRI.
      Magnetic Resonance in Med. 2013 Jun;69(6):1534–1540.
    This is the main function, which will invoke the optimization and
    produce the directions incrementally.
    '''
    Ks  = np.array(Ks)
    K   = np.sum(Ks)
    S   = len(Ks)
    Ksc = np.zeros(S, dtype=int) # current directions per shell, will eventually match Ks
    vectors = [] # spherical coordinates, only theta and phi
    shells  = []
    W  = np.full((S,S), 1-alpha)
    np.fill_diagonal(W, alpha)
    for k in range(K):
        shell   = next_shell(Ksc, Ks)
        vector  = find_one_vector(shell, vectors, shells, W, maxiter=maxiter)
        vectors.append(vector)
        shells.append(shell)
        Ksc[shell] += 1
    vectors = np.array(vectors) # spherical coordinates, theta and phi only
    r       = np.ones((len(vectors),1))
    vectors = np.hstack((r, vectors)) # now with the radii
    vectors = sph2cart(vectors) # cartesian coordinates
    shells  = np.array(shells)[:,None] # shells start counting from 0
    return vectors, shells

def find_one_vector(shell, vectors, shells, W, maxiter=500, epsilon=1e-9): # --
    '''
    Find an optimal new direction, given the existing ones in their
    corresponding shells. Since we are doing incrementally, we don't need to
    optimize all vectors simulaneously, and can replace SLS-QP for L-BFGS-B,
    dispensing with equality constraints by working with spherical coordinates,
    and dispensing also with the need for the full normalizing weights from
    Eqns 2 and 4 of the paper (W depends only on alpha).  Makes it simpler
    and, as it turns out, also faster.
    Vectors are in spherical coordinates (theta and phi only; no radius), so
    that the bounds are well defined (0,pi) and (0,2*pi).
    '''
    init_guess = cart2sph(random_vectors(1))[0,1:] # theta and phi only, as r = 1
    result     = scopt.minimize(
        cost_one_vector,
        init_guess,
        args    = (shell, vectors, shells, W),
        method  = 'L-BFGS-B',
        bounds  = [(0,np.pi), (0,2*np.pi)],
        options =  {'maxiter': maxiter,
                    'ftol'   : epsilon,
                    'disp'   : False})
    vector = result.x # theta and phi only
    return vector

def cost_one_vector(vector, shell, vectors, shells, W, epsilon=1e-9): # ----
    '''
    Cost function for adding a new vector to an existing set.
    In this function, inputs vector and vectors are in spherical coordinates
    '''
    u = sph2cart(np.array((1, *vector)))
    Energy = 0.0
    for i, v in enumerate(vectors):
        v       = sph2cart(np.array((1, *v)))
        diff    = np.sum((u - v) ** 2) + epsilon
        ssq     = np.sum((u + v) ** 2) + epsilon
        energy  = 1/diff + 1/ssq  # Eqn 3 of the paper
        w       = W[shell, shells[i]]
        Energy += w * energy
    return Energy

def next_shell(Ksc, Ks): # ----------------------------------------------------
    '''
    Define the next shell to work on based on which shell has currently
    a greater "deficit" in number of directions. The deficit is given by the
    max deviation from the desired proportion of directions across shells
    Ks  : Desired number of directions per shell
    Ksc : Current number of directions per shell
    This task could be done faster by taking some operations to outside the
    loop, but it would make it all less clear, so let's leave them here
    '''
    K  = np.sum(Ks)  # desired total number of directions
    Kc = np.sum(Ksc) # current total number of directions
    Ps = Ks / K      # desired proportion of directions across shells
    if Kc == 0:
        deficits = Ps
    else:
        Pc       = Ksc / Kc # current proportion
        deficits = Ps - Pc
    nextshell = np.argmax(deficits) # starts at 0
    return nextshell


# ######[ MISCELLANEOUS FUNCTIONS ]############################################

def directions_per_shell(K, S, distrib='linearly'): # =========================
    '''
    Given a number of directions and shells, distribute the direction counts
    across shells so that they are evenly distributed (same count across all 
    shells), linearly distributed (linearly more as we move to outer shells),
    or quadratically distributed (quadratically more as we move to outer shells).
    '''
    if   distrib in [0, 'uniform', 'uniformly', 'evenly', 'constant', 'constantly']:
        Ps = np.ones(S)
    elif distrib in [1, 'lin', 'linear', 'linearly']:
        Ps = np.arange(1,S+1, dtype=float)
    elif distrib in [2, 'quad', 'quadratic', 'quadratically']:
        Ps = np.arange(1,S+1, dtype=float) ** 2
    Ks0 = Ps/np.sum(Ps)*K
    Ks  = Ks0.astype(int)
    while np.sum(Ks) < K:
        Kdiff = Ks0 - Ks
        Ks[np.argmax(Kdiff)] += 1
    while np.sum(Ks) > K:
        Kdiff = Ks0 - Ks
        Ks[np.argmin(Kdiff)] -= 1
    return Ks

def random_vectors(K): # ======================================================
    '''
    Generate K random vectors uniformly distributed on the
    surface of the unit sphere.
    '''
    vectors = np.random.randn(K,3)
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:,None]
    return vectors

def cart2sph(cart): # =========================================================
    '''
    Converts cartesian (x,y,z) coordinates to spherical (r,theta,phi) coords.
    Theta is measured in relation +z, and phi in relation to +x
    Azimuth = phi, Elevation = pi/2 - theta
    '''
    r     = np.linalg.norm(cart, axis=1)
    theta = np.arccos(cart[:,2] / r)
    theta = np.where(r == 0, 0, theta) # rare case in which r = 0
    phi   = np.mod(np.arctan2(cart[:,1], cart[:,0]), 2*np.pi)
    sph   = np.column_stack((r, theta, phi))
    return sph

def sph2cart(sph): # ==========================================================
    '''
    Converts spherical (r,theta,phi) coordinates to cartesian (x,y,z) coords.
    Theta is measured in relation +z, and phi in relation to +x
    Azimuth = phi, Elevation = pi/2 - theta
    '''
    r, theta, phi = sph.T
    cart = np.column_stack((
        r * np.sin(theta) * np.cos(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(theta) ))
    return cart

def intersperse_rows(A, B): # =================================================
    '''
    Generic function to intersperse rows of matrices A and B evenly.
    The first row of the merged matrix will always come from the smallest.
    If they have the same size, the first row comes from B.
    '''
    rA, cA = A.shape
    rB, cB = B.shape
    rM     = rA + rB
    if cA != cB:
        raise ValueError('The two matrices must have the number of columns.')
    M = np.empty((rM,cA))
    if rA == rB:
        M[0::2] = B
        M[1::2] = A
        return M
    elif rA > rB:
        nI   = int(np.ceil(rM/rB)*rB)
        idx  = np.reshape(np.arange(nI), newshape=(rB,int(nI/rB)))
        idxA = np.ravel(idx[:,1:], order='C')[:rA]
        idxB = idx[:,0]
    elif rA < rB:
        nI   = int(np.ceil(rM/rA)*rA)
        idx  = np.reshape(np.arange(nI), newshape=(rA,int(nI/rA)))
        idxA = idx[:,0]
        idxB = np.ravel(idx[:,1:], order='C')[:rB]
    M[idxA,:] = A
    M[idxB,:] = B
    return M

def add_b0(vectors, bvalues, K0=1, where='interspersed', unique=False): # =====
    '''
    Add a number K0 of bvalues=0 to a pair of vectors and bvalues.
    Some vendors require directions to be unique (e.g., Philips), so we'll
    create random vectors (that are unlikely to be repeated) for bvalues=0.
    '''
    if K0 == 0:
        return vectors, bvalues
    if unique:
        vectors0 = random_vectors(K0)
    else:
        vectors0 = np.zeros((K0, 3))
    bvalues0 = np.zeros((K0, 1))
    if where == 'start':
        vectors = np.concatenate((vectors0, vectors), axis=0)
        bvalues = np.concatenate((bvalues0, bvalues), axis=0)
    elif where == 'end':
        vectors = np.concatenate((vectors, vectors0), axis=0)
        bvalues = np.concatenate((bvalues, bvalues0), axis=0)
    elif where == 'interspersed':
        vectors = intersperse_rows(vectors, vectors0)
        bvalues = intersperse_rows(bvalues, bvalues0)
    return vectors, bvalues

def remove_b0(bvalues, vectors, shells): # ====================================
    '''
    Remove vectors and shells that have a corresponding bvalue = 0.
    Return also an index that allows putting these back.
    '''
    idx     = bvalues != 0
    bvalues = bvalues[idx]
    vectors = vectors[idx]
    shells  = shells[idx]
    return bvalues, vectors, shells, idx

def scale_vectors(vectors, bvalues): #=========================================
    '''
    Scale the norm of the vectors for scanners that don't take bvalues for
    each direction, just a maximum bvalue (e.g., Siemens, GE).
    '''
    bmax    = bvalues.max()
    if bmax == 0:
        vectors = np.full(vectors.shape,0)
    else:
        vectors = np.sqrt(bvalues/bmax)*vectors
    return vectors

def unscale_vectors(vectors, bmax): #==========================================
    '''
    Scake back to unit norm the vectors from scanners that don't take bvalues 
    for each direction, just a maximum bvalue (e.g., Siemens, GE), and use the
    max bvalue to figure out what the bvalues should be.
    '''
    norm     = np.linalg.norm(vectors, axis=1)
    bvalues  = np.round(norm**2 * bmax).astype(int)[:,None]
    idx = norm != 0
    vectors[idx,:] /= norm[idx,None]
    return vectors, bvalues

# ######[ PLOTTING FUNCTIONS ]#################################################

def geodesic_sphere(n=2): # ===================================================
    '''
    Make a geodesic sphere based on an icosahedron, with recursive subdivision.
    '''
    def subdivide(vtx, fac):
        vtxn = vtx.tolist()
        facn = []
        edg2mid = {}
        next_index = len(vtx)
        def get_midpoint(i, j):
            nonlocal next_index
            if i > j:
                i, j = j, i
            key = (i,j)
            if key not in edg2mid:
                v_i, v_j = vtx[i], vtx[j]
                mid = (v_i + v_j) / 2
                mid /= np.linalg.norm(mid)
                vtxn.append(mid)
                edg2mid[key] = next_index
                next_index += 1
            return edg2mid[key]
        for f in fac:
            a, b, c = f
            m_ab = get_midpoint(a, b)
            m_bc = get_midpoint(b, c)
            m_ca = get_midpoint(c, a)
            facn.extend([
                [a, m_ab, m_ca],
                [m_ab, b, m_bc],
                [m_ca, m_bc, c],
                [m_ab, m_bc, m_ca]])
        return np.array(vtxn), np.array(facn)
    
    g = (1 + np.sqrt(5)) / 2 # Golden ratio
    vtx = np.array([
        [0, 1, g], [0, -1, g], [0, 1, -g], [0, -1, -g],
        [1, g, 0], [-1, g, 0], [1, -g, 0], [-1, -g, 0],
        [g, 0, 1], [g, 0, -1], [-g, 0, 1], [-g, 0, -1]
    ])
    norm = np.linalg.norm(vtx, axis=1, keepdims=True)
    vtx /= norm
    fac = np.array([
            [ 0,  1,  8], [ 0, 10,  1], [ 2,  9,  3], [ 2,  3, 11],
            [ 1,  7,  6], [ 3,  6,  7], [ 5,  0,  4], [ 2,  5,  4],
            [ 6,  9,  8], [ 8,  9,  4], [ 7, 10, 11], [ 5, 11, 10],
            [ 1,  6,  8], [ 0,  8,  4], [10,  7,  1], [10,  0,  5],
            [ 6,  3,  9], [ 7, 11,  3], [ 2,  4,  9], [ 2, 11,  5] ])
   
    for _ in range(n):
        vtx, fac = subdivide(vtx, fac)
    return vtx, fac

def uv_sphere(nlat=20, nlon=20): # ============================================
    '''
    Make a UV (latitude-longitude) sphere.
    '''
    vtx = []
    vtx.append([0.0, 0.0, 1.0]) # North pole
    for i in range(1, nlat): # Rings between poles
        theta = i * np.pi / nlat
        for j in range(nlon):
            phi = j * 2 * np.pi / nlon
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            vtx.append([x, y, z])
    vtx.append([0.0, 0.0, -1.0]) # South pole
    
    # Faces (all triangular)
    fac = []
    for j in range(nlon): # North pole
        v0 = 0  # North pole
        v1 = 1 + j
        v2 = 1 + (j + 1) % nlon
        fac.append([v0, v1, v2])
    for k in range(1, nlat - 1): # Middle section
        for j in range(nlon):
            v0 = 1 + (k - 1) * nlon + j
            v1 = 1 + (k - 1) * nlon + (j + 1) % nlon
            v2 = 1 + k * nlon + (j + 1) % nlon
            v3 = 1 + k * nlon + j
            # Split quad into two triangles
            fac.append([v0, v1, v2])
            fac.append([v0, v2, v3])
    N = 1 + (nlat - 1) * nlon  # South pole index
    for j in range(nlon): # South pole
        v0 = N
        v1 = 1 + (nlat - 2) * nlon + (j + 1) % nlon
        v2 = 1 + (nlat - 2) * nlon + j
        fac.append([v0, v1, v2])
    vtx = np.array(vtx)
    fac = np.array(fac)
    return vtx, fac

def northern_hemisphere(vtx, fac): # ==========================================
    '''
    Select vertices and faces of the northern hemisphere (z >= 0) from a sphere.
    '''
    # Vertex indices of interest
    zpos   = vtx[:,2] >= 0
    newvtx = vtx[zpos]
    
    # Map from old to new indices
    zidx   = np.where(zpos)[0]
    newidx = np.full(len(vtx), -1, dtype=int)
    newidx[zidx] = np.arange(len(zidx))
    
    # Select faces where all vertices have z >= 0
    fpos   = np.all(zpos[fac], axis=1)
    fidx   = fac[fpos]
    newfac = newidx[fidx]
    return newvtx, newfac

def plot_directions(vectors, shells=None, bvalues=None, 
                    filename=None, colorby='shell', sphere='uv',
                    style='quiver', reproject=False): # =======================
    """
    Plot a multishell with hemispheres for each shell or or qspace direction
    scheme. 'style' can be 'quiver' for arrows, 'points' for scatter, 'blobs'
    for little circles representing each direction, 'rings' for the circles
    without filling, or 'qspace' for a scatter in the cartesian coordinate
    system.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    
    # Do not open a window if saving to a file
    if filename is not None:
        matplotlib.use('agg')  # Force Agg backend, override existing
    
    # Here we want the shell indices to start at 1
    if shells is None and bvalues is not None:
        if not reproject:
            vectors = scale_vectors(vectors, bvalues)
        _, shells = np.unique(bvalues, return_inverse=True)
        shells    = shells + 1
        bmax      = bvalues.max()
        bvalues1  = bvalues / bmax # scaled to between 0 and 1
    elif shells is not None and bvalues is None:
        shells    = shells - shells.min() + 1
        bmax      = None
        bvalues1  = shells / shells.max() # fake b-values, and scaled to between 0 and 1
        if not reproject:
            vectors = vectors * shells/shells.max()
    else:
        raise ValueError('Must provide either "shells" or "bvalues", not both')
    uB1 = np.unique(bvalues1)
    uS  = np.unique(shells)
    S   = len(uS)
    K   = len(vectors)

    # We only need one hemisphere of the shells; let's flip so that they
    # are all on the +z hemisphere, i.e., above the "equator".
    idx = vectors[:,2] < 0
    vectors[idx] = -vectors[idx]
    
    # Custom colormap . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if 'multishell' not in plt.colormaps():
        colors = [(.95, .05, .05),
                  (.05, .85, .05),
                  (.05, .05, .95)]
        cmap = mcolors.LinearSegmentedColormap.from_list('multishell', colors, N=2**10)
        matplotlib.colormaps.register(cmap)
    cmap   = plt.get_cmap('multishell')
    if   colorby == 'shell':
        C     = S
        Cidx  = shells
        uCidx = uS
    elif colorby == 'acquisition':
        C     = K
        Cidx  = np.arange(K) + 1
        uCidx = Cidx
    colors = [cmap(c / (C - 1)) if C > 1 else cmap(0.5) for c in range(C)]
    cmap   = mcolors.ListedColormap(colors)
    bounds = np.concatenate([[uCidx[0]   - 0.5], 
                             (uCidx[:-1] + uCidx[1:]) / 2, 
                             [uCidx[-1]  + 0.5]])
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)
    
    # Set up axes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection='3d')
    if style == 'qspace':
        ax.grid(True)
        ax.xaxis.pane.set_alpha(.5)
        ax.yaxis.pane.set_alpha(.5)
        ax.zaxis.pane.set_alpha(.5)
        ax.xaxis.line.set_visible(True)
        ax.yaxis.line.set_visible(True)
        ax.zaxis.line.set_visible(True)
    else:
        ax.grid(False)
        ax.xaxis.pane.set_alpha(0)
        ax.yaxis.pane.set_alpha(0)
        ax.zaxis.pane.set_alpha(0)
        ax.xaxis.line.set_visible(False)
        ax.yaxis.line.set_visible(False)
        ax.zaxis.line.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(0,1) 
    ax.set_box_aspect([1,1,.5])
    ax.set_proj_type('ortho')
    
    # Directions  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if style == 'points':
        sc = ax.scatter(vectors[:,0], vectors[:,1], vectors[:,2], 
                        c=Cidx.ravel(order='C'), cmap=cmap, norm=norm, alpha=1, s=5)
    elif style == 'qspace':
        sphere = None
        sc = ax.scatter(vectors[:,0], vectors[:,1], vectors[:,2], 
                        c=Cidx.ravel(order='C'), cmap=cmap, norm=norm, alpha=1, s=5)
    elif style == 'quiver':
        for i, vec in enumerate(vectors):
            if reproject:
                arrow_length_ratio = .1
            else:
                if bvalues is None:
                    arrow_length_ratio = .1/bvalues1[i]
                else:
                    arrow_length_ratio = .1/np.sqrt(bvalues1[i])
            ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], 
                      color=cmap(norm(Cidx[i])), alpha=1, linewidth=1.5,
                      arrow_length_ratio=arrow_length_ratio)
    elif style == 'rings':
        for i, vec in enumerate(vectors):
            if reproject:
                r = 0.05
            else:
                if bvalues is None:
                    r = 0.05 * bvalues1[i]
                else:
                    r = 0.05 * np.sqrt(bvalues1[i])
            ref   = np.zeros(3, dtype=float)
            ref[np.argmin(np.abs(vec))] = 1.0
            u     = np.cross(vec, ref)
            u    /= np.linalg.norm(vec)
            v     = np.cross(vec, u)
            v    /= np.linalg.norm(vec)
            phi   = np.linspace(0, 2*np.pi, 36)
            circle = vec + r * (np.cos(phi)[:,None]*u + np.sin(phi)[:,None]*v)
            ax.plot(circle[:,0], circle[:,1], circle[:,2], color=cmap(norm(Cidx[i])))
    elif style == 'blobs':
        from mpl_toolkits.mplot3d import art3d
        for i, vec in enumerate(vectors):
            if reproject:
                r = 0.05
            else:
                if bvalues is None:
                    r = 0.05 * bvalues1[i]
                else:
                    r = 0.05 * np.sqrt(bvalues1[i])
            ref   = np.zeros(3, dtype=float)
            ref[np.argmin(np.abs(vec))] = 1.0
            u     = np.cross(vec, ref)
            u    /= np.linalg.norm(vec)
            v     = np.cross(vec, u)
            v    /= np.linalg.norm(vec)
            phi   = np.linspace(0, 2*np.pi, 36)
            circle = vec + r * (np.cos(phi)[:,None]*u + np.sin(phi)[:,None]*v)
            ax.add_collection(art3d.Poly3DCollection(circle[None,:], facecolors=cmap(norm(Cidx[i])), linewidth=0))
        ax.view_init(elev=90, azim=0)
    else:
        raise ValueError("Style must be 'points', 'quiver', 'rings', 'blobs', or 'qspace'.")
    
    # Reference sphere  . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if sphere is not None:
        if sphere in ['geo', 'geodesic']:
            vtx, fac = geodesic_sphere(n=2)
        elif sphere in ['uv', 'latlong']:
            vtx, fac = uv_sphere(nlat=40, nlon=80) # choose an even number for nlat, otherwise the north hemi will be incomplete near the equator
            vtx, fac = northern_hemisphere(vtx,fac)
            edg = set()
            for f in fac:
                edg.update([(min(a,b), max(a,b)) for a, b in zip(f, np.roll(f,-1))])
            for b1 in uB1:
                if reproject:
                    svtx = vtx
                else:
                    if bvalues is None:
                        svtx = vtx * b1
                    else:
                        svtx = vtx * np.sqrt(b1)
                tri = [svtx[f,:] for f in fac]
                poly = Poly3DCollection(tri, facecolors='gray', edgecolors='none', alpha=0.1)
                ax.add_collection3d(poly)
                x, y, z = [], [], []
                for i, j in edg:
                    x.extend([svtx[i,0], svtx[j,0], np.nan])
                    y.extend([svtx[i,1], svtx[j,1], np.nan])
                    z.extend([svtx[i,2], svtx[j,2], np.nan])
                ax.plot(x, y, z, color='gray', alpha=0.1, linewidth=0)
                if reproject:
                    break

    # Colorbar  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if   colorby == 'shell':
        label = 'Shell number'
        ticks = uCidx
        ticklabels = [f'{c:.0f}' for c in uCidx]
    elif colorby == 'acquisition':
        label = 'Acquisition order'
        ticks = np.linspace(1,K,3)
        ticklabels = [f'{c:.0f}' for c in (1,(1+K)/2,K)]
    if style == 'points':
        plt.colorbar(sc, ax=ax, label=label, ticks=ticks, fraction=0.05, shrink=0.5)
    else:
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(uCidx)
        cbar = plt.colorbar(sm, ax=ax, label=label, ticks=ticks, fraction=0.05, shrink=0.5)
        cbar.set_ticklabels(ticklabels)
    ax.set_clip_on(True)
    
    # Save figure . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
        plt.close(fig)


# ######[ PARSE ARGUMENTS ]####################################################

def take_args(args, defaults): # ==============================================
    '''
    Parse arguments
    '''
    # Arguments and options will be stored in the 'state' dict
    state = {
        'new'        : None,      # create a new set?
        'method'     : None,      # simultaneous or incremental
        'K'          : None,      # total number of directions
        'S'          : None,      # number of shells
        'distr'      : None,      # how directions are distributed across shells
        'Ks'         : None,      # list with directions per shell (alternative to supplying K, S, and distr)
        'reorder'    : None,      # whether do the optimal reordering
        'addb0'      : None,      # list with 3 integers
        'removeb0'   : None,      # remove existing b0?
        'bvalues'    : None,      # list with b-values
        'bmax'       : None,      # max b-value, needed to figure out the bvalues for some scanner formats
        'in'         : None,      # input file to be read
        'intype'     : None,      # format of the input file
        'out'        : None,      # output file to be saved
        'metadata'   : None,      # metadata for Siemens scanners
        'comment'    : None,      # comment for Philips scanners
        'vectors'    : None,      # vectors that will be saved/plotted
        'shells'     : None }     # shells that will be saved/plotted
    state = state | defaults

    # Loop over arguments
    a = 0
    while a < len(args):
        
        if args[a] == '--help':
            # Help text
            a += 1
            
        elif args[a] == '--new':
            # Information to create a new set of directions distributed across shells
            state['method'] = args[a+1]
            if args[a+2][0] in '[(':
                state['Ks'] = eval(args[a+2])
                state['K']  = sum(state['Ks'])
                state['S']  = len(state['Ks'])
                a += 3
            else:
                state['K']     = eval(args[a+2])
                state['S']     = eval(args[a+3])
                state['distr'] = args[a+4]
                a += 5
        
        elif args[a] == '--alpha':
            # Alpha parameter that weights intra shell (alpha=1) vs
            # inter-shell (alpha=0) energy
            state['alpha'] = float(eval(args[a+1]))
            a += 2
        
        elif args[a] == '--reorder':
            # Do the optimal reordering for simulaneous or file conversion?
            a += 1
        
        elif args[a] == '--bvalues':
            # Bvalues, as needed to convert between certain file formats
            state['bvalues'] = eval(args[a+1])
            a += 2
            
        elif args[a] == '--bmax':
            # Bvalues, as needed to convert between certain file formats
            state['bmax'] = eval(args[a+1])
            a += 2
            
        elif args[a] == '--addb0':
            # Add bvalues = 0?
            state['addb0'] = eval(args[a+1])
            a += 2
            
        elif args[a] == '--removeb0':
            # Do the optimal reordering for simulaneous or file conversion?
            a += 1
            
        elif args[a] == '--load':
            # File that will be read with information
            state['in']     = args[a+1]
            state['intype'] = args[a+2].lower()
            a += 3
        
        elif args[a] == '--save':
            # File that will be saved with results
            state['outprefix'] = args[a+1]
            if len(args) > a+2 and args[a+2][0] != '-':
                state['outtype'] = args[a+2]
                a += 3
            else:
                a += 2

        elif args[a] == '--plot':
            # Output figures
            if len(args) > a+1 and args[a+1][0] != '-':
                state['plotformat'] = args[a+1].lower()
                a += 2
            else:
                a += 1
        else:
            raise ValueError('Unknown option: {}'.format(args[a]))
    return args, state

def check_errors(args, state): # ==============================================
    '''
    Check for most likely errors when calling the command.
    '''
    
    # These output formats include or not bvalues, explicitly or implicitly through rescaling
    inclbvals = ('fsl', 'mrtrix3', 'siemens', 'ge', 'philips')
    nobvals   = ('caruyer')
    
    # These formats scale the vectors and don't explicitly include bvals and thus need bmax
    needbmax  = ('siemens', 'ge')
    
    if '--new' in args and '--load' in args:
        raise ValueError('Cannot use --new and --load in the same execution.')

    if '--new' in args:
        
        if state['method'] not in ('simultaneous', 'incremental'):
            raise ValueError('Unknon method {}. Must be "simultaneous" or "incremental".'.format(state['method']))
            
        if state['Ks'] is None and (state['K'] is None or state['S'] is None or state['distr'] is None):
            raise ValueError(
                'To create new shells, must specify either:'
                ' (a) a list with directions for each shell or '
                ' (b) number of directions, number of shells, '
                '     and how directions should be distributed.')
            
        if state['K'] is not None and state['S'] is not None:
            if state['S'] > state['K']:
                raise ValueError('Cannot have more shells than directions.')
            
        if state['bvalues'] is not None and state['S'] is not None:
            if (len(state['bvalues']) != state['S']):
                raise ValueError('Number of bvalues and number of shells must match.')
    
    if '--alpha' in args:
        
        if ('--new' not in args) and ('--reorder' not in args):
            print('Option --alpha has no effect without --new or --reorder. Ignoring it.')
        elif state['alpha'] < 0 or state['alpha'] > 1:
            raise ValueError('The value provided with --alpha must be betwee 0 and 1.')
    
    if '--load' in args:
        
        if state['intype'] is None:
            raise ValueError('Must specify an input file format with --load.')
            
        if state['intype'] in needbmax:
            if state['bmax'] is None:
                raise ValueError('Must specify bmax to load a {} file format.'.format(state['intype']))
        else:
            if state['bmax'] is not None:
                raise ValueError('Cannot provide --bmax with {} input format.'.format(state['intype']))
                
        if state['intype'] in ('fsl', 'mrtrix3', 'siemens', 'ge', 'philips'):
            if state['bvalues'] is not None:
                raise ValueError('Cannot provide --bvalues with {} input format.'.format(state['intype']))
    
    if '--save' in args:
        if (state['outtype'] in inclbvals) or ('--addb0' in args) or ('--removeb0' in args):
            if  ('--new'  in args) or \
                ('--load' in args and state['intype'] in nobvals):
                    if '--bvalues' not in args:
                        raise ValueError('Unknown bvalues. Provide them with --bvalues.\n'
                                         'Example: --bvalues "[1000,2000,3000]"')
    
    if '--plot' in args:
        if state['plotformat'] not in ('interactive', 'pdf', 'png', 'jpeg', 'webp', 'svg'):
            raise ValueError('Unknown option to plot.\n'
                             'Examples: --plot interactive\n'
                             '          --plot pdf')
    return args, state


# ######[ MAIN, GLORIOUS FUNCTION ]############################################
if __name__ == "__main__":
    
    # Defaults  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    defaults = {
            'alpha'      : 0.75,           # unlike in the paper, we use alpha=0.75 (not 0.50)
            'maxiter'    : 500,            # maximum number of iterations for the optimizers
            'outtype'    : 'caruyer',      # format of the output file
            'plotformat' : 'interactive' } # format to save outputs (pdf or png)
    
    # Parse arguments and check for errors
    args, state = take_args(sys.argv[1:], defaults)
    args, state = check_errors(args, state)
    
    # Call help . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if '--help' in args or len(args) == 0:
        print_help()
        sys.exit(0)
        
    # Create a new optimal set . . . . . . . . . . . . . . . . . . . . . . . . 
    if '--new' in args:
        
        print('Creating direction set with {} optimization'.format(state['method']))
        if not state['Ks']:
            state['Ks'] = directions_per_shell(state['K'], state['S'], distrib=state['distr'])
        if state['method'] == 'simultaneous':
            state['vectors'], state['shells'] = \
                simultaneous_optimization(state['Ks'], alpha=state['alpha'], maxiter=state['maxiter'])
        elif state['method'] == 'incremental':
            state['vectors'], state['shells'] = \
                incremental_optimization( state['Ks'], alpha=state['alpha'], maxiter=state['maxiter'])
        
    # Load an existing set  . . . . . . . . . . . . . . . . . . . . . . . . . . 
    if '--load' in args:
        
        print('Reading input file: {}'.format(state['in']))
        if   state['intype'] == 'fsl'  or state['intype'] == 'bids' or \
             state['intype'] == 'bvec' or state['intype'] == 'bval':
            # Read from FSL format
            state['vectors'], state['bvalues'] = read_fsl(state['in'])
            _, shells          = np.unique(state['bvalues'], return_inverse=True)
            state['shells']    = shells[:,None]
            
        elif state['intype'] == 'mrtrix3':
            # Read from MRtrix3 format
            state['vectors'], state['bvalues'] = read_mrtrix3(state['in'])
            _, shells          = np.unique(state['bvalues'], return_inverse=True)
            state['shells']    = shells[:,None]
            
        elif state['intype'] == 'siemens':
            # Read directions from Siemens scanners
            state['vectors'], state['metadata'] = read_siemens(state['in'])
            state['vectors'], state['bvalues']  = unscale_vectors(state['vectors'], state['bmax'])
            _, shells          = np.unique(state['bvalues'], return_inverse=True)
            state['shells']    = shells[:,None]
            
        elif state['intype'] == 'ge':
            # Read directions from GE scanners
            state['vectors']   = read_ge(state['in'])
            state['vectors'], state['bvalues'] = unscale_vectors(state['vectors'], state['bmax'])
            _, shells          = np.unique(state['bvalues'], return_inverse=True)
            state['shells']    = shells[:,None]
        
        elif state['intype'] == 'philips':
            # Read directions from Philips scanners
            state['vectors'], state['bvalues'], state['comment'] = read_philips(state['in'])
            _, shells          = np.unique(state['bvalues'], return_inverse=True)
            state['shells']    = shells[:,None]
            
        elif state['intype'] == 'caruyer':
            # Read directions from Caruyer's webtool format
            state['vectors'], state['shells'] = read_caruyer(state['in'])
    
    # Ensure bvalues is a column vector, count number of shells . . . . . . . .
    if state['bvalues'] is not None:
        if type(state['bvalues']) is list:
            state['bvalues']   = np.array(state['bvalues'])
            state['bvalues']   = state['bvalues'][state['shells'][:,0]][:,None]
    state['S'] = len(np.unique(state['shells'][:,0]))
    
    # Optimally reorder the directions if asked . . . . . . . . . . . . . . . .
    if '--reorder' in args:
        
        if state['method'] == 'incremental':
            print('Directions already in optimal order. Skipping option "--reorder".')
        else:
            print('Reordering directions')
            if state['bvalues'] is None:
                bidx = slice(None)
            else:
                bidx = state['bvalues'][:,0] != 0
            state['vectors'][bidx], state['shells'][bidx], idx = \
                optimal_reordering(state['vectors'][bidx], state['shells'][bidx])
            if state['bvalues'] is not None:
                state['bvalues'][bidx] = state['bvalues'][bidx][idx]

    # Remove definitively bvalues=0 if the user asked . . . . . . . . . . . . .
    if '--removeb0' in args:
        
        print('Removing existing bvalues = 0')
        bidx             = state['bvalues'][:,0] != 0
        state['vectors'] = state['vectors'][bidx]
        state['shells']  = state['shells'][bidx]
        state['bvalues'] = state['bvalues'][bidx]
    
    # Add bvalues = 0 if the user requested . . . . . . . . . . . . . . . . . .
    if '--addb0' in args:
        
        print('Adding bvalues = 0')
        
        # Keep the existing shells for b>0, reindex from 0 if needed
        bidx    = state['bvalues'][:,0] != 0
        shellsb = state['shells'][bidx]
        shellsb = shellsb - shellsb.min()
        
        # Add new bvalues = 0
        if state['outtype'] == 'philips':
            unique = True
        else:
            unique = False
        state['vectors'], state['bvalues'] = \
            add_b0(state['vectors'], state['bvalues'], \
                   K0=state['addb0'][1], unique=unique, where='interspersed')
        state['vectors'], state['bvalues'] = \
            add_b0(state['vectors'], state['bvalues'], \
                   K0=state['addb0'][0], unique=unique, where='start')
        state['vectors'], state['bvalues'] = \
            add_b0(state['vectors'], state['bvalues'], \
                   K0=state['addb0'][2], unique=unique, where='end')
        
        # Restore shells
        bidx = state['bvalues'][:,0] != 0
        state['shells'] = np.zeros(state['bvalues'].shape)
        state['shells'][bidx] = shellsb + 1
    
    # Save direction set . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    if '--save' in args:
        
        print('Saving directions in file with prefix: {}'.format(state['outprefix']))
        if   state['outtype'] == 'fsl'  or state['outtype'] == 'bids' or \
             state['outtype'] == 'bvec' or state['outtype'] == 'bval':
            # Write in FSL format                
            write_fsl(state['outprefix'],
                      state['vectors'],
                      state['bvalues'])
            
        elif state['outtype'] == 'mrtrix3':
            # Write in MRtrix3 format
            write_mrtrix3(state['outprefix'] + '.b',
                          state['vectors'],
                          state['bvalues'])
            
        elif state['outtype'] == 'siemens':
            # Write directions for Siemens scanners
            vectors = scale_vectors(state['vectors'], state['bvalues'])
            write_siemens(state['outprefix'] + '.dvs',
                          vectors,
                          metadata=state['metadata'])
            
        elif state['outtype'] == 'ge':
            # Write directions for GE scanners
            vectors = scale_vectors(state['vectors'], state['bvalues'])
            write_ge(state['outprefix'] + '.dat',
                     vectors)
        
        elif state['outtype'] == 'philips':
            # Write directions for Philips scanners
            if state['bvalues'][0,0] != 0:
                print('For Philips, the first bvalue must be 0. Adding it now.')
                rndvec = random_vectors(1)
                state['vectors'] = np.vstack((rndvec,state['vectors']))
                state['bvalues'] = np.vstack((0,state['bvalues']))
            b0idx = state['bvalues'][:,0] == 0
            if np.unique(state['vectors'][b0idx], axis=0).shape[0] != state['vectors'][b0idx].shape[0]:
                print('For Philips, there cannot be repeated directions, not even for b0. ' + 
                      'Replacing repeated directions for b0 with random directions.')
                state['vectors'][b0idx] = random_vectors(sum(b0idx))
            write_philips(state['outprefix'] + '.txt',
                          state['vectors'],
                          state['bvalues'],
                          comment=state['comment'])
            
        elif state['outtype'] == 'caruyer':
            # Write directions in Caruyer's webtool format
            if state['bvalues'] is None:
                bidx = np.full((len(state['vectors']),), True)
            else:
                bidx = state['bvalues'][:,0] != 0
            write_caruyer(state['outprefix'] + '.txt',
                          state['vectors'][bidx],
                          state['shells'][bidx])

    # Plot . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    if '--plot' in args:

        # Remove bvalues=0 as we won't plot these
        if state['bvalues'] is None:
            bidx = slice(None)
        else:
            bidx = state['bvalues'][:,0] != 0
        
        # No need to reproject if there is only 1 shell, as the figures are identical
        if state['S'] == 1:
            Reproject = (True,)
        else:
            Reproject = (False, True)

        # Make the figures for the interactive or save mode
        # In the interactive mode, the user saves through the plot window
        # In the non-interactive mode, the images are saved as is
        # If bvalues are available, directions and shells are scaled according
        # to the bvalues. Otherwise, they are scaled proportionally to 1, 2, 3, etc.
        if state['plotformat'] == 'interactive':
            print('Plotting directions in interactive mode')
            for style in ['quiver', 'blobs']:
                for reproject in Reproject:
                    for colorby in ('shell','acquisition'):
                        if state['bvalues'] is None:
                            plot_directions(state['vectors'][bidx],
                                            shells    = state['shells'][bidx],
                                            style     = style,
                                            colorby   = colorby,
                                            reproject = reproject)
                        else:
                             plot_directions(state['vectors'][bidx],
                                             bvalues   = state['bvalues'][bidx],
                                             style     = style,
                                             colorby   = colorby,
                                             reproject = reproject)   
        else:
            print('Plotting directions and saving image filesas: {}-*.{}'.format(state['outprefix'], state['plotformat']))
            for style in ['quiver', 'blobs']:
                for reproject in Reproject:
                    for colorby in ('shell','acquisition'):
                        filename = '{}-{}-reproj{}-colorby{}.{}'.format(state['outprefix'], style, reproject, colorby, state['plotformat'])
                        if state['bvalues'] is None:
                            plot_directions(state['vectors'][bidx],
                                            shells    = state['shells'][bidx],
                                            filename  = filename,
                                            colorby   = colorby,
                                            style     = style,
                                            reproject = reproject)
                        else:
                            plot_directions(state['vectors'][bidx],
                                            bvalues   = state['bvalues'][bidx],
                                            filename  = filename,
                                            colorby   = colorby,
                                            style     = style,
                                            reproject = reproject)
    print('Done!')
    sys.exit(0)