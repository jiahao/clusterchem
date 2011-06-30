#!/usr/bin/env python

from twisted.internet import protocol, reactor
from twisted.internet.defer import Deferred
from twisted.internet.utils import getProcessOutput

import os

#Lower-level I/O handler
class MyHandler(protocol.ProcessProtocol):
    def connectionMade(self):
        print 'ohai', self.transport

    def childDataReceived(self, childFD, chunk):
        print 'Some Sub Data:', childFD, repr(chunk)

    def processEnded(self, reason):
        print 'ended because:', str(reason)

#Higher-level job handler
class SGEWaiter:
    """This class will submit a program to the Sun Grid Engine,
    then wait for it to complete."""

    def __init__(self, myreactor = None):
        """@param reactor a Twisted reactor"""
        self.reactor = myreactor
        if myreactor is None: #Use global reactor
            self.reactor = reactor
        
        self.handler = MyHandler()
        self.MainLoop()

    def stage1(self):
        cmd = '/bin/bash'
        print 's1!'
        #self.reactor.spawnProcess(self.handler, cmd, args = [cmd, '-c', """
        output = getProcessOutput(cmd, args = ['-c', """
echo Stage1
sleep 3
echo Stage1 done
"""], env = os.environ, reactor = self.reactor)
        print 's1!'
        output.addCallbacks(self.stage2)
        return output

    def stage2(self, s1result):
        print 's2!'
        cmd = '/bin/sh'
        output = getProcessOutput(cmd, args = ['-c', """
echo Stage2
sleep 3
echo Stage2 done
"""], env = os.environ, reactor = self.reactor)
        print 's2!'
        print s1result
        output.addCallbacks(self.tearDown)
        return output

    def setUp(self):
        print 'Setting up'

    def tearDown(self, s2result):
        """Clean up after itself."""
        print s2result
        print 'Tearing down'
        #self.reactor.stop()

    def MainLoop(self):
        self.setUp()
        d1 = self.stage1()
        #def stage1ready(stage1output):
        #    d2 = self.stage2()
        #    def stage2ready(stage2output):
        #        # stage1output and stage2output in scope!
        #        self.tearDown()
        #    d2.addCallback(stage2ready)
        #    return d2
        #
        #d1.addCallback(stage1ready)
        #return d1

if __name__ == '__main__':
    SGEWaiter()

    #Force terminate after 10s
    reactor.callLater(10, reactor.stop)
    reactor.run()

