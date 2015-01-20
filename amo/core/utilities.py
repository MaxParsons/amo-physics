'''
Created on Jan 19, 2015

@author: Max
'''
from datetime import datetime


class simulations(object):
    @staticmethod
    def talk(instr):
        now = datetime.now()
        print '[TALK] ' + now.strftime("%H:%M:%S.%f").rstrip('0') + '--- ' + instr
