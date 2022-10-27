#!/usr/bin/env python
# coding: utf-8

import os
import sys

file = sys.argv[1]

basename = os.path.basename(file)
month = basename[-11:-4]
print(month)
