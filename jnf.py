from sympy import *

def compute_J_and_X(A, lamda):

    # get the dimension n of the nxn matrix A
    n = A.shape[0]
    JList = []
    XList = []
    for eigenV in lamda:
        N = A - (eigenV * eye(n))
        #creating a basis
        B = N.nullspace()
        #dimension of the basis
        currentLen = len(B)
        oldLen = -1
        M = N
        # in dimensions we save 
        # how many NEW dimensions we have at each level when multiplying by N so
        # it also tells us how many chains we have that reach that level
        dimensions = [currentLen]
        kernels = [B]
        while (currentLen != oldLen):
            M *= N
            oldLen = currentLen
            currentLen = len(M.nullspace())
            dimensions.append(currentLen - oldLen)
            kernels.append(M.nullspace())
        # eliminate the 0 since there is no new dimension at the end 
        # (the difference between last level and the previous is just zero)
        dimensions.pop()
        # create J
        # in jordanBoxSize we save mj = dj - dj+1 (how many chains end in each level)
        # so the index + 1 tells us what size the jordanBox has to be + 1 size 
        # the sum of the values in jordanBoxSize is the geometric multiplicity 
        jordanBoxSize = []
        for i in range(len(dimensions)):
            if i < len(dimensions) - 1:
                jordanBoxSize.append(dimensions[i] - dimensions[i + 1])
            else:
                jordanBoxSize.append(dimensions[i])
        for i in reversed(range(len(jordanBoxSize))):
            # get j boxes of size i+1
            for j in range(jordanBoxSize[i]):
                jordanBlock = Matrix.jordan_block(i + 1, eigenV)
                JList.append(jordanBlock)
        # create X
        # dictionary to get accepted vectors for each level (basis)
        acceptedVectorsPerLevel = {}
        for i in range(len(kernels)):
            acceptedVectorsPerLevel[i] = []
        for level in reversed(range(len(jordanBoxSize))):
            # we get the number of chains that reach that level so we know
            # how often we have to iterate over that
            numberOfChains = jordanBoxSize[level]
            if level > 0:
                previousBasis = kernels[level - 1]
            else:
                previousBasis = []
            # iterate k times so if we have 2 chains here as an example
            # we have to get 2 startingvectors for that level 
            for k in range(numberOfChains):
                # create referenceBasis with all vectors that already exist on
                # previous level so we do not get linearly dependent vectors
                referenceBasis = previousBasis + acceptedVectorsPerLevel[level]
                # we check each vector on our current level
                for vector in kernels[level]:
                    # just take the first vector if the referenceBasis is empty 
                    if(len(referenceBasis) == 0):
                        startVector = vector
                        acceptedVectorsPerLevel[level].append(startVector)
                        break
                    # compare referenceBasis and current vector to check with
                    # basis that is one level below
                    TestMatrix = Matrix.hstack(*(referenceBasis + [vector]))
                    referenceMatrix = Matrix.hstack(*referenceBasis)
                    if (TestMatrix.rank() == referenceMatrix.rank() + 1):
                        startVector = vector
                        acceptedVectorsPerLevel[level].append(startVector)
                        break
                # create matrix with n * (level + 1) that saves our startVector
                # of our current chain on the right and then calculate
                # our new vectors with x(level - 1) = N * x(level) 
                length = level + 1
                matrix = zeros(n, length)
                matrix[:, -1] = startVector
                for i in reversed(range(length - 1)):
                # calculate new vector from the column next to it by N
                    matrix[:, i] = N * matrix[:, i+1]
                    acceptedVectorsPerLevel[i].append(matrix[:, i])
                XList.append(matrix)

    J = diag(*JList)
    X = Matrix.hstack(*XList)
    return J, X
