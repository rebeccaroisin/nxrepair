#!/usr/bin/env python

# Copyright (c) 2014, rebeccaroisin
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import math
import random


class IntervalNode( object ):
    def __init__( self, start, end, linenum=0, other=None ):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.linenum = linenum
        self.other = other
    def insert( self, start, end, linenum=0, other=None ):
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert( start, end, linenum, other )
            else:
                self.right = IntervalNode(start, end, linenum, other )
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert( start, end, linenum, other )
            else:
                self.left = IntervalNode(start, end, linenum, other )
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left: 
            root.maxend = max( root.end, root.right.maxend, root.left.maxend )
            root.minend = min( root.end, root.right.minend, root.left.minend )
        elif root.right: 
            root.maxend = max( root.end, root.right.maxend )
            root.minend = min( root.end, root.right.minend )
        elif root.left:
            root.maxend = max( root.end, root.left.maxend )
            root.minend = min( root.end, root.left.minend )
        return root

    def rotateright( self ):
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root
        
    def rotateleft( self ):
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root

    def intersect( self, start, end, report_func ):
        if start < self.end and end > self.start: report_func( self )
        if self.left and start < self.left.maxend:
            self.left.intersect( start, end, report_func )
        if self.right and end > self.start:
            self.right.intersect( start, end, report_func )

    def traverse( self, func ):
        if self.left: self.left.traverse( func )
        func( self )
        if self.right: self.right.traverse( func )