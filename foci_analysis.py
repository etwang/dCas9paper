import os, subprocess, operator
import matplotlib
matplotlib.use('Agg')
from pylab import *
from skimage import img_as_float, morphology, filters, measure, io, img_as_uint, img_as_ubyte
import time
import Mysbatch
params = {'legend.fontsize': 8}
rcParams.update(params)



# Open a zeiss image file.
def openFile(img_f):
    from czifile import CziFile
    with CziFile(img_f) as czi:
        image_arrays = czi.asarray()
    print 'Shape:', image_arrays.shape
    # Returns several dimensions:
    # ?, colors, ?, Z, X, Y, ?
    return image_arrays


# Save a part of an image as a numpy array file.
def sliceImage(img_f, xmin, xmax, ymin, ymax, out_f):
   
    xmin, xmax, ymin, ymax = map(int, [xmin, xmax, ymin, ymax])
 
    image_arrays = openFile(img_f)
    images = image_arrays[0, :, 0, :, :, :, 0]
    colors, zsize, xsize, ysize = images.shape
    save(out_f, images[:, :, xmin:xmax, ymin:ymax])


def selectSubsetsWrapper(imagedir, xdim, ydim, sliceddir, clustertype):

    files = [f for f in os.listdir(imagedir) if f.endswith(".czi")] 
    for f in files:
        outname = os.path.join(sliceddir, f.split(".")[0])
        cmd = 'python %s/script.py dCas9_v3.selectSubsets %s %s %s %s'%(\
            sys.path[0],\
            os.path.join(imagedir, f), xdim, ydim, outname)
        print cmd
        if clustertype == 'slurm':
            scriptOptions = {'ppn':1, 'jobname': 'selectSubsets',\
                'memory': '40gb', 'qos':'ewang-b', 'time':'2:00:00'} 
            Mysbatch.launchJob(cmd, scriptOptions, verbose=True) 



# Save only the parts of images that have a fiber.
def selectSubsets(img_f, xdim, ydim, outname):
   
    from scipy import stats 
    xdim, ydim = map(int, [xdim, ydim]) 

    # Open file
    print 'Opening image'
    image_arrays = openFile(img_f)
    images = image_arrays[0, :, 0, :, :, :, 0]
    colors, zsize, xsize, ysize = images.shape

    # Record dapi and gfp arrays
    print 'Computing threshold'
    dapi = images[0, :, :, :]
    gfp = images[1, :, :, :]

    # Use GFP to threshold the fiber
    thresh = filters.threshold_otsu(gfp)
    print thresh
     
    xcoords = [x for x in range(xsize) if x%xdim == 0]
    xcoords.append(xsize)
    ycoords = [y for y in range(ysize) if y%ydim == 0]
    ycoords.append(ysize)

    print 'Iterating through image'
    for x in range(len(xcoords) - 1):
        for y in range(len(ycoords) - 1):
            print xcoords[x], ycoords[y]
            block = gfp[:, xcoords[x] : xcoords[x + 1],\
                ycoords[y] : ycoords[y + 1]]
            mask = block > thresh
            mask = morphology.binary_opening(mask, ones((5, 5, 5)))
            try:
                zmin, zmax, xmin, xmax, ymin, ymax = bbox2_3D(mask)
            except:
                zmin, zmax, xmin, xmax, ymin, ymax = [0, 0, 0, 0, 0, 0]
            if (zmax - zmin) * (xmax - xmin) * (ymax - ymin) >\
                0.05 * zsize * xdim * ydim:
                print 'Use me'
                xmin = xmin + xcoords[x]
                xmax = xmax + xcoords[x]
                ymin = ymin + ycoords[y]
                ymax = ymax + ycoords[y]
                dapi_f = outname + ".z%s-%s.x%s-%s.y%s-%s.dapi.npy"%(\
                    zmin, zmax, xmin, xmax, ymin, ymax)
                save(dapi_f, dapi[zmin : zmax, xmin : xmax, ymin : ymax])
                gfp_f = outname + ".z%s-%s.x%s-%s.y%s-%s.gfp.npy"%(\
                    zmin, zmax, xmin, xmax, ymin, ymax)
                save(gfp_f, gfp[zmin : zmax, xmin : xmax, ymin : ymax])


def showSubsets(indir, prefix, out_f):
    
    files = [f for f in os.listdir(indir) if f.startswith(prefix) and \
        'dapi' in f]
    xvals = []
    yvals = []
    for f in files:
        print f
        pre, z, x, y = f.split(".")[:4]
        xmin, xmax = map(int, x[1:].split("-"))
        ymin, ymax = map(int, y[1:].split("-"))
        xvals.append(xmax)    
        yvals.append(ymax)    
    im = zeros((max(xvals), max(yvals)))
    
    for f in files:
        print f
        pre, z, x, y = f.split(".")[:4]
        xmin, xmax = map(int, x[1:].split("-"))
        ymin, ymax = map(int, y[1:].split("-"))
        
        image = load(os.path.join(indir, f))
        im[xmin : xmax, ymin : ymax] = image.max(axis=0)

    imshow(im)
    colorbar()
    axis('off')
    savefig(out_f)


def bbox2_3D(img):

    r = any(img, axis=(1, 2))
    c = any(img, axis=(0, 2))
    z = any(img, axis=(0, 1))

    rmin, rmax = where(r)[0][[0, -1]]
    cmin, cmax = where(c)[0][[0, -1]]
    zmin, zmax = where(z)[0][[0, -1]]

    return rmin, rmax, cmin, cmax, zmin, zmax

def segmentAndProcessWrapper(indir, outdir, clustertype):
    
    files = [f for f in os.listdir(indir)]
    for f in files:
        out_f = os.path.join(outdir, f + ".png")
        cmd = 'python %s/script.py dCas9_v3.segmentAndProcess %s %s'%(\
            sys.path[0], os.path.join(indir, f), out_f)
        print cmd
        if clustertype == 'slurm':
            scriptOptions = {'ppn':1, 'jobname': 'segment',\
                'memory': '5gb', 'qos':'ewang-b', 'time':'1:00:00'} 
            Mysbatch.launchJob(cmd, scriptOptions, verbose=True) 


def segmentAndProcess(image_f, out_f):
    from scipy import stats 

    print 'Opening image'
    if image_f.endswith(".czi"):
        image_arrays = openFile(image_f)
        images = image_arrays[0, :, 0, :, :, :, 0]
    else:
        images = load(image_f)
    print 'Done opening image'
   
    dapi = images[0, :, :, :]
    dapi = img_as_float(dapi)
   
    figure(figsize=(8.5, 11))
 
    # Show raw DAPI
    subplot2grid((3, 2), (0, 0))
    title('Raw DAPI')
    imshow(dapi.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off') 

    dapi = filters.gaussian(dapi, 3)    
    subplot2grid((3, 2), (1, 0))
    title('Gaussian Filter')
    imshow(dapi.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off') 

    subplot2grid((3, 2), (2, 0))
    title('Histogram')
    hist(dapi.flatten(), color='#CCCCCC')
    m = dapi.mean()
    sd = dapi.std()
    axvline(m, color='k')
    axvline(m + 5 * sd, linestyle='--', color='k')
    xlim(0, 1)
    ylim(0, 1e5)

    binary = zeros(dapi.shape)
    binary[where(dapi > m + 5 * sd)] = 1
    subplot2grid((3, 2), (0, 1))
    title('Thresholded')
    imshow(binary.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off')

    from scipy import ndimage as ndi

    print 'Cleaning'
    binary = ndi.binary_fill_holes(binary)
    binary = morphology.binary_opening(binary, morphology.ball(2))
    #binary = morphology.binary_closing(binary, morphology.ball(4))
    subplot2grid((3, 2), (1, 1))
    title('Cleaned')
    imshow(binary.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off')

    labels = measure.label(binary)
    labels = morphology.remove_small_objects(labels, 1000)

    for region in measure.regionprops(labels):
        if region.area > 20000:
            labels[where(labels == region.label)] = 0
        

    subplot2grid((3, 2), (2, 1))
    title('Cleaned & Labeled')
    imshow(labels.max(axis=0), cmap='Spectral')
    colorbar()
    axis('off')
    
    suptitle(image_f)
    savefig(out_f)

def segmentAndProcessWithGFPWrapper(indir, outdir, clustertype):
    
    files = [f for f in os.listdir(indir)]
    for f in files:
        out_f = os.path.join(outdir, f + ".png")
        txt_f = os.path.join(outdir, f + ".txt")
        cmd = 'python %s/script.py dCas9_v3.segmentAndProcessWithGFP %s %s %s'%(\
            sys.path[0], os.path.join(indir, f), out_f, txt_f)
        print cmd
        if clustertype == 'slurm':
            scriptOptions = {'ppn':1, 'jobname': 'segment',\
                'memory': '5gb', 'qos':'ewang-b', 'time':'1:00:00'} 
            Mysbatch.launchJob(cmd, scriptOptions, verbose=True) 


def segmentAndProcessWithGFP(image_f, out_f, txt_f):
    from scipy import stats 

    print 'Opening image'
    if image_f.endswith(".czi"):
        image_arrays = openFile(image_f)
        images = image_arrays[0, :, 0, :, :, :, 0]
    else:
        images = load(image_f)
   
    dapi = images[0, :, :, :]
    dapi = img_as_float(dapi)
    gfp = images[1, :, :, :]
    gfp = img_as_float(gfp)
   
    figure(figsize=(11, 8.5))
 
    # Show raw DAPI
    subplot2grid((3, 4), (0, 0))
    title('Raw DAPI')
    imshow(dapi.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off') 

    # Gaussian filter
    print 'Gaussian Filter'
    dapi = filters.gaussian(dapi, 3)    
    subplot2grid((3, 4), (1, 0))
    title('Gaussian Filter')
    imshow(dapi.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off') 

    # Histogram for DAPI
    subplot2grid((3, 4), (2, 0))
    hist(dapi.flatten(), color='#CCCCCC', normed=True)
    m = dapi.mean()
    sd = dapi.std()
    axvline(m, color='k')
    axvline(m + 5 * sd, linestyle='--', color='k')
    xlim(0, 1)
    ylim(0, .1)
    xlabel('DAPI intensity')
    ylabel('Frequency')

    print 'Thresholding'
    binary = zeros(dapi.shape)
    binary[where(dapi > m + 5 * sd)] = 1
    subplot2grid((3, 4), (0, 1))
    title('Thresholded')
    imshow(binary.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off')

    from scipy import ndimage as ndi

    print 'Cleaning'
    binary = ndi.binary_fill_holes(binary)
    binary = morphology.binary_opening(binary, morphology.ball(2))
    #binary = morphology.binary_closing(binary, morphology.ball(4))
    subplot2grid((3, 4), (1, 1))
    title('Cleaned')
    imshow(binary.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off')

    print 'Labeling'
    labels = measure.label(binary)
    labels = morphology.remove_small_objects(labels, 1000)

    numobj = 0
    for region in measure.regionprops(labels):
        if region.area > 20000:
            labels[where(labels == region.label)] = 0
        else:
            numobj += 1
        
    subplot2grid((3, 4), (2, 1))
    title('Further cleaned &\nLabeled')
    imshow(labels.max(axis=0), cmap='Spectral')
    xlabel('%s objects'%(numobj))
    colorbar()
    axis('off')
    
    # Show raw GFP 
    subplot2grid((3, 4), (0, 2))
    title('Raw GFP')
    imshow(gfp.max(axis=0), cmap='Greys_r', vmin=0, vmax=1) 
    colorbar()
    axis('off') 


    print 'Getting GFP background'
    # Get GFP background
    subplot2grid((3, 4), (1, 2))
    gfpbg = 3 * morphology.opening(gfp, morphology.ball(3))
    title('GFP Opened')
    imshow(gfpbg.max(axis=0), cmap='Greys_r', vmin=0, vmax=1) 
    colorbar()
    axis('off') 
    
    print 'Subtracting GFP background'
    sub = gfp - gfpbg
    sub[where(sub < 0)] = 0
    # Get GFP background
    subplot2grid((3, 4), (2, 2))
    title('GFP subtracted')
    imshow(sub.max(axis=0), cmap='Greys_r', vmin=0, vmax=1) 
    colorbar()
    axis('off') 

    print 'Foci in Nuclei'
    foci = zeros(gfp.shape)
    foci[where(labels > 0)] = sub[where(labels > 0)]
    subplot2grid((3, 4), (0, 3))
    title('Foci in nuclei')
    imshow(foci.max(axis=0), cmap='Greys_r') 
    colorbar()
    axis('off') 

    # Histogram for GFP 
    subplot2grid((3, 4), (1, 3))
    hist(sub.flatten(), 20, color='k', normed=True, histtype='step', label='All')
    hist(foci[where(labels > 0)].flatten(), 20, color='r', normed=True, histtype='step', label='In nuclei')
    xlim(0, 1)
    ylim(0, .1)
    legend(loc='upper right')
    xlabel('GFP intensity')
    ylabel('Frequency')

    print 'Quantitating Intensities'
    # Get intensity per nucleus
    intensities = []
    txt = open(txt_f, 'w')
    txt.write("#Nucleus\tArea\tIntensity\tPerPixel\n")
    for region in measure.regionprops(labels):
        idx = where(labels == region.label)
        perpixel = foci[idx].sum() / region.area
        intensities.append(perpixel)
        txt.write("\t".join(map(str, [region.label, region.area, foci[idx].sum(), perpixel])) + "\n")
    txt.close()
    subplot2grid((3, 4), (2, 3))
    hist(intensities, 20, color='#CCCCCC')
    xlabel('Int. per pixel')
    ylabel('Frequency')
    xticks(rotation=90)
    xlim(0, .1)

    suptitle(image_f)
    tight_layout()
    subplots_adjust(top=.9)

    savefig(out_f, dpi=300)


# Plot max intensity stack from array.
def maxIntensity(imagedir, xdim, ydim, channel, outdir):

    channel, xdim, ydim = map(int, [channel, xdim, ydim]) 
    
    for f in os.listdir(imagedir):
        if f.endswith(".czi"):
            out_f = os.path.join(outdir, f.split(".")[0] + ".png")
            npy_f = os.path.join(outdir, f.split(".")[0] + ".npy")
            if not os.path.exists(out_f):

                image_arrays = openFile(os.path.join(imagedir, f))
                images = image_arrays[0, :, 0, :, :, :, 0]
                colors, zsize, xsize, ysize = images.shape

                print 'Plotting image', f
                dapi = images[channel, :, :, :].max(axis=0)

                figure()
                imshow(dapi.T, cmap='Greys', origin='lower')

                xcoords = [x for x in range(xsize) if x%xdim == 0]
                xcoords.append(xsize)
                ycoords = [y for y in range(ysize) if y%ydim == 0]
                ycoords.append(ysize)
                for x in xcoords:
                    axvline(x)
                for y in ycoords:
                    axhline(y)
                ylim(0, ysize)
                xlim(0, xsize)
                xticks(rotation=90)
                colorbar()
                        
                print 'Saving'
                savefig(out_f)
                save(npy_f, dapi)

def plotSlices(slice_f, thumbdir, out_f):
   
    sampleToCoords = {} 
    for line in open(slice_f):
        if not line.startswith("#"):
            sample, x, y = line.strip().split()
            x = float(x)
            y = float(y)
            if sample not in sampleToCoords:
                sampleToCoords[sample] = []
            sampleToCoords[sample].append([x, y])

    import matplotlib.patches as patches
    n = 1
    figure(figsize=(8.5, 11))
    samples = sampleToCoords.keys()
    samples.sort()
    for sample in samples:
        print sample
        ax = subplot(6, 3, n) 
        img = load(os.path.join(thumbdir, sample + ".npy"))
        imshow(img.T, cmap='Greys', origin='lower')
        data = sampleToCoords[sample]
        for x, y in data:
            ax.add_patch(patches.Rectangle((x, y), 500, 500, alpha=.1))
        title(sample)
        xlim(0, img.shape[0]) 
        ylim(0, img.shape[1]) 
        axis('off')
        n += 1
    savefig(out_f)

def saveSlices(slice_f, imagedir, outdir):        
   
    sampleToCoords = {} 
    for line in open(slice_f):
        if not line.startswith("#"):
            sample, x, y = line.strip().split()
            x = int(x)
            y = int(y)
            if sample not in sampleToCoords:
                sampleToCoords[sample] = []
            sampleToCoords[sample].append([x, y])

    samples = sampleToCoords.keys()
    samples.sort()
    for sample in samples:
        print sample
        data = sampleToCoords[sample]
        image_arrays = openFile(os.path.join(imagedir, sample + ".czi"))
        img = image_arrays[0, :, 0, :, :, :, 0]
        for x, y in data:
            block = img[:, :, x : x + 500, y : y + 500]
            out_f = os.path.join(outdir, sample + ".x%s.y%s.npy"%(x, y))
            save(out_f, block)

def matchDistributions(indir, samplelist, out_f):
   
    from stats import ks_2samp
    import seaborn as sns

    samples = []
    for line in open(samplelist):
        sample, obj = line.strip().split()
        samples.append(sample)

    sampleToData = {}
    for f in os.listdir(indir):
        sample = f.split(".")[0]
        if sample in samples and f.endswith(".txt"):
            for line in open(os.path.join(indir, f)):
                if not line.startswith("#"):
                    nuc, area, intensity, pixel = line.strip().split("\t")
                    try:
                        sampleToData[sample].append(float(pixel))
                    except:
                        sampleToData[sample] = [float(pixel)]
       
    figure(figsize=(40, 40)) 
    sns.set_style("white")
    
    from scipy import stats 
    samples = sampleToData.keys()
    samples.sort()
    for i in range(len(samples)):
        print i
        for j in range(i + 1, len(samples)):
            y1, x1, p1 = hist(sampleToData[samples[i]], linspace(0, .1, 20))
            y2, x2, p2 = hist(sampleToData[samples[j]], linspace(0, .1, 20))
            subplot2grid((len(samples), len(samples)), (j, i))
            scatter(y1, y2)
            xlabel(samples[i])
            ylabel(samples[j])
    savefig(out_f)


def plotDistributions(indir, samplelist, out_f):
   
    from stats import ks_2samp
    import seaborn as sns

    ntype = {'Uninf':[], 'empty':[], 'CAG':[]} 
    samples = []
    for line in open(samplelist):
        samples.append(line.strip().split()[0])

    for f in os.listdir(indir):
        ntypename = f.split("_")[0]
        if f.split(".")[0] in samples and f.endswith(".txt"):
            for line in open(os.path.join(indir, f)):
                if not line.startswith("#"):
                    nuc, area, intensity, pixel = line.strip().split("\t")
                    ntype[ntypename].append(float(pixel))
       
    figure(figsize=(2.5, 1.8)) 
    sns.set_style("white")
   
    nTypeToLabel = {'Uninf': 'Uninfected',\
        'empty': 'dSaCas9$_{ctrl}$',\
        'CAG': 'dSaCas9$_{CAG}$'}
 
    from scipy import stats 
    for nuc in ['Uninf', 'empty', 'CAG']:
        #y, x, p = hist(ntype[nuc], 20,\
        #    visible=False, normed=True)
        #plot(x[1:], y, '-', label=nuc + ' (%s nuclei)'%(len(ntype[nuc])))
        vmin = stats.scoreatpercentile(ntype[nuc], 10)
        vmax = stats.scoreatpercentile(ntype[nuc], 90)
        print nuc, vmin, vmax 
        sns.kdeplot(array(ntype[nuc]), label=nTypeToLabel[nuc] + \
            ' (%s nuclei)'%(len(ntype[nuc])))

    
    d, p = ks_2samp(ntype['Uninf'], ntype['CAG'])
    print 'Uninf', 'CAG', d, p
    d, p = ks_2samp(ntype['empty'], ntype['CAG'])
    print 'Empty', 'CAG', d, p
    d, p = ks_2samp(ntype['empty'], ntype['Uninf'])
    print 'Empty', 'Uninf', d, p

    yticks(linspace(0, 80, 5), fontsize=6)
    xticks(linspace(0, .075, 4), linspace(0, 3, 4), fontsize=6)
    ylim(0, 80)
    xlim(0, .075)
  
    params0 = {'legend.fontsize': 6}
    rcParams.update(params0) 
    legend(loc='upper right')
    xlabel('Mean Intensity per unit volume\nin nucleus (a. u.)', fontsize=8)
    ylabel('Frequency', fontsize=8)
    subplots_adjust(left=.2, bottom=.3)
    sns.despine()
    savefig(out_f)


# Subtract background of fibers and do histogram equalization for visualization
def cleanForFigure(czi_f, coords_f, dapi_f, gfp_f):

    coords = []  
    for line in open(coords_f):
        fname, xmin, xmax, ymin, ymax = line.strip().split()
        coords.append(map(int, [xmin, xmax, ymin, ymax]))

 
    image_arrays = openFile(czi_f)
    images = image_arrays[0, :, 0, :, :, :, 0]
    colors, zsize, xsize, ysize = images.shape
    dapi = images[0, :, :, :]
    gfp = images[1, :, :, :]

    newdapi = zeros((zsize, xsize, ysize), dtype=float)
    newgfp = zeros((zsize, xsize, ysize), dtype=float)
   
    for xmin, xmax, ymin, ymax in coords: 
        print xmin, xmax, ymin, ymax

        # Subtract backgrounds.
        dapiblock = img_as_float(dapi[:, xmin : xmax, ymin : ymax])

        gfpblock = img_as_float(gfp[:, xmin : xmax, ymin : ymax])
        
        print 'GFP Opening'
        gfpbg = 3 * morphology.opening(gfpblock, morphology.ball(3))
        
        newdapi[:, xmin : xmax, ymin : ymax] = dapiblock 
        newgfp[:, xmin : xmax, ymin : ymax] = gfpblock - gfpbg
    
    print dapi.shape
    io.imsave(dapi_f, img_as_ubyte(newdapi.max(axis=0)))
    io.imsave(gfp_f, img_as_ubyte(newgfp.max(axis=0)))




