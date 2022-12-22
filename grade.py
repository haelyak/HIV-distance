import unittest
import dist as hw


class SimpleTests(unittest.TestCase):

    def testB(self):
        # propDifferent
        self.assertTrue(hw.propDifferent('A-T','CTT')==0.5,msg="problem in propDifferent")
        
    def testC(self):
        # propDifferent
        self.assertTrue(round(hw.propDifferent('GGG','TTG'),4)==0.6667,msg="problem in propDifferent")

    def testD(self):
        # propDifferent
        self.assertTrue(hw.propDifferent('CGA','C-A')==0.0,msg="problem in propDifferent")

    def testE(self):
        # propDifferent
        self.assertTrue(round(hw.propDifferent('------ACCGTACGATCAGTA','ACGTACGATGCAG--------'),4)==0.7143,msg="problem in propDifferent")

    def testF(self):
        # distances
        D=hw.distances(['a','b'],['A','A'])
        self.assertTrue(type(D)==type({}),msg="problem in distances")
        self.assertTrue(len(D.keys())==4,msg="problem in distances")
        self.assertTrue(type(list(D.keys())[0])==type(()),msg="problem in distances")
        self.assertTrue(type(list(D.values())[0])==type(0),msg="problem in distances")

    def testG(self):
        # distances
        D=hw.distances(['a','b','c'],['ACTT','A-GT','TATT'])
        self.assertTrue(round(D[('a','b')],4)==0.4408,msg="problem in distances")
        
    def testH(self):
        # distances
        D=hw.distances(['a','b','c'],['ACTT','A-GT','TATT'])
        self.assertTrue(round(D[('c','b')],4)==1.6479,msg="problem in distances")

    def testI(self):
        # distances
        D=hw.distances(['a','b','c'],['ACTT','A-GT','TATT'])
        sortedL=[(('a', 'a'), 0), (('a', 'b'), 0.4408399986765892), (('a', 'c'), 0.8239592165010822), (('b', 'a'), 0.4408399986765892), (('b', 'b'), 0), (('b', 'c'), 1.6479184330021643), (('c', 'a'), 0.8239592165010822), (('c', 'b'), 1.6479184330021643), (('c', 'c'), 0)]
        self.assertTrue(compareDistD(D,sortedL),msg="problem in distances")
        
    def testJ(self):
        # distances
        D=hw.distances(['w','x','y','z'],['ACGGAC','ACGTG-','-CGGGA','ACGGTA'])
        L=sorted(D.items())

def compareDistD(D,sortedL):
    L=sorted(D.items())
    if len(L) != len(sortedL): return False
    else:
        for i in range(len(L)):
            if L[i][0]!=sortedL[i][0] or round(L[i][1],4)!=round(sortedL[i][1],4):
                return False
        return True
        
if __name__ == '__main__':
    unittest.main()
