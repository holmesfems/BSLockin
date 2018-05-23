import sys, time

class progbar:
    def __init__(self, period=100, bars=32):
        self._period  = period
        self.bars     = bars
        self.active   = True
        self.start    = time.time()

    def dispose(self):
        if self.active:
            self.active = False
            self.update(self._period)
            sys.stdout.write("\n")

    def __del__(self):
        self.dispose()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.dispose()

    def period(self):
        return self._period

    def update(self, tick):
        rate = tick / self._period

        # progress rate
        str = "{0:4d}% ".format(int(rate*100))

        # progress bar
        bar_prog = int(rate * self.bars)
        str += "|"
        str += "#" * (            bar_prog)
        str += "-" * (self.bars - bar_prog)
        str += "|"

        # calc end
        elapsed = time.time() - self.start
        predict = (elapsed * (1 - rate) / rate) if not rate == 0.0 else 0
        s = int(predict)
        m, s = divmod(s, 60)
        h, m = divmod(m, 60)
        #d, h = divmod(h, 24)
        str += "{0:3d}:{1:02d}:{2:02d}\r".format(h,m,s)

        sys.stdout.write(str)
        sys.stdout.flush()

def main():
    # how to use
    with progbar(500) as pb:
        for i in range(pb.period()):
            time.sleep(0.01)
            if i % 10 == 0:
                pb.update(i)

if __name__ == '__main__':
    main()

