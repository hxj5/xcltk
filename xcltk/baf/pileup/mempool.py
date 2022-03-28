# Memory Pool
# Author: Xianjie Huang

# TODO: implement push_back().
class MemPool:
    """Memory Pool
    @note Input should be either a @class unit, which implements reset() method, or
          two functions @func func_init and @func func_reset.
    """
    def __init__(self, unit = None, func_init = None, func_reset = None):
        self.unit = unit
        self.func_init = func_init
        self.func_reset = func_reset
        self.data = []
        self.l = 0
        self.n = 0
        self.is_reset = 0

    def destroy(self):
        self.data = []

    def get(self):
        if self.l >= self.n:
            if self.n >= len(self.data):
                n = 16 if len(self.data) == 0 else len(self.data)
                self.data.extend([None for _ in range(n)])
            self.data[self.n] = self.unit() if self.unit else self.func_init()
            self.n += 1
            self.l += 1
            return(self.data[self.l - 1])
        else:
            if self.is_reset:
                if self.unit:
                    self.data[self.l].reset()
                else:
                    self.func_reset(self.data[self.l])
            self.l += 1
            return(self.data[self.l - 1])

    def reset(self):
        self.l = 0
        self.is_reset = 1

class PoolTestClass:
    def __init__(self):
        self.a = 10

    def reset(self):
        self.a = 10

def __test_init():
    return({"key":20})

def __test_reset(x):
    x["key"] = 20

if __name__ == "__main__":
    print("test func input ...")
    pool = MemPool(None, __test_init, __test_reset)
    x = pool.get()
    print("x is %d" % x["key"])
    print("pool: %s" % str(pool.data))
    x["key"] -= 1
    print("with modification, x is %d" % x["key"])
    print("pool: %s" % str(pool.data))
    print("reset pool ...")
    pool.reset()
    x = pool.get()
    print("x is %d" % x["key"])
    print("pool: %s" % str(pool.data))
    x["key"] -= 1
    print("with modification, x is %d" % x["key"])
    print("pool: %s" % str(pool.data))

    print()
    
    print("test class input ...")
    pool = MemPool(PoolTestClass)
    x = pool.get()
    print("x is %d" % x.a)
    print("pool: %s" % str([d.a if d else None for d in pool.data]))
    x.a -= 1
    print("with modification, x is %d" % x.a)
    print("pool: %s" % str([d.a if d else None for d in pool.data]))
    print("reset pool ...")
    pool.reset()
    x = pool.get()
    print("x is %d" % x.a)
    print("pool: %s" % str([d.a if d else None for d in pool.data]))
    x.a -= 1
    print("with modification, x is %d" % x.a)
    print("pool: %s" % str([d.a if d else None for d in pool.data]))
    
