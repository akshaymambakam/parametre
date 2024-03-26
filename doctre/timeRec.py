import time, os, sys
start_time = time.time()
os.system(sys.argv[1])
print("--- Total time taken: %s seconds ---" % (time.time() - start_time))
