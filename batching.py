# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
import os, time, sys
import multiprocessing as mp
import subprocess

def worker(py, fn):
    try:
        job = subprocess.Popen(['python', py, fn])     
        job.wait()
    except OSError as e:
        print e, py, fn
    return fn

def callback(fn):
    print fn, 'done'
    
def apply_async_with_callback(py, filelist):
    cpu = mp.cpu_count()-2
    
    print '# cpus: %d' % cpu       
    pool = mp.Pool(processes = cpu)
    for f in filelist:
        pool.apply_async(worker, (py, f), callback=callback)
    pool.close()
    pool.join()
    
if __name__ == '__main__':
    """
    Batch process laz/las files in a specific folder with specific python script
    usage: $ python batching.py <script.py> <folder location>
    example: $ python batching.py dtm.py ../example_data/laz/ (! note the final slash)
    """    
    py = 'biomass.py'
    path = '../../example_data/laz/'
    
    #py = sys.argv[1]
    #path = sys.argv[2]
    
    filelist = os.listdir(path) 
    fn =  []
    t0 = time.time()
    for f in filelist: 
        if (f[-3:] == 'las') or (f[-3:] == 'laz'):    
            fn.append(path+f)
    apply_async_with_callback(py, fn)
