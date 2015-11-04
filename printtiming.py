import time

def print_timing(func):
    def wrapper(*arg, **kwarg):
        t1 = time.time()
        localtime = time.asctime( time.localtime(time.time()) )
        print "\n  :::::::::: %s start time: %s ::::::::::\n" % (func.func_name, localtime)
        res = func(*arg, **kwarg)
        t2 = time.time()
        localtime2 = time.asctime( time.localtime(time.time()) )
        print "\n  :::::::::: %s end time: %s ::::::::::" % (func.func_name, localtime2)
        print '  :::::::::: %s took %0.2f seconds (%.2f min).::::::::::\n' % (func.func_name, (t2-t1)*1.0, (t2-t1)/60.0)
        return res
    return wrapper