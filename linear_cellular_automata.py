"""
linear_cellular_automata.py
author: moatjon

This file contains functions to compute and display a linear cellular automaton (LCA).
"""

import numpy as np


def printGeneration(aGeneration):
    myprintString = ""
    for i in aGeneration:
        if i == True:
            myprintString += "*"
        else:
            myprintString += " "
    print(myprintString)


def getSumMatrix(n):
    """
    Generates summation matrix for universal linear cellular automaton usage
    
    arguments:
    n -- the height and width of the summation matrix

    return:
    sumMat -- a summation matrix of shape n x n

    Product of the summation matrix and a generation vector is a vector with sums of
    neighbourhood values, where the left neighbour has a value of four, the center
    cell has a value of two and the right neighbour has a value of one. Thus the result
    can easily be matched against a lookup table which describes automaton behavior.
    """

    # Create vector for the three diagonals which are to be filled with values 1, 2 and 4
    oneVec = [1] * (n - 1)
    twoVec = [2] * n
    fourVec = [4] * (n - 1)
    # Create the matrix by summation of the vectors at different diagonals (k)
    sumMat = np.diag(fourVec, k=-1) + np.diag(twoVec) + np.diag(oneVec, k=1)
    return sumMat


def lutify(indices, lut):
    """
    Parses a list of lookup table (LUT) indexes according to a given lookup table

    arguments:
    indices -- a list of indices that correspond to certain values in a LUT
    lut -- a LUT against which indices are checked for output values

    returns:
    lutified -- a list of values corresponding to the passed indices
    """
    
    lutified = []
    # For every value in indices, find which value corresponds in the LUT
    for index in indices:
        lutified.append(lut[index])
    return lutified

        
def propagateLCA(predecessor, nGenerations, automatonLut):
    """
    Calculates propagated LCA generations based on an automaton truth table and summation matrix

    arguments:
    predecessor -- a first generation which serves as the basis for the next
    nGenerations -- the number of generations to calculate, excluding predecessor
    automatonLut -- a lookup table (LUT) which describes automaton behavior

    returns:
    generations -- a list of calculated generations, starting with the predecessor
    """
    
    generations = [predecessor]
    generationLength = len(predecessor)
    # Generate a summation matrix for easy calculation of LUT indices
    sumMatrix = getSumMatrix(generationLength)
    for i in range(nGenerations):
        # Calculate the next generation based on the most recently calculated one
        nextGen = lutify(sumMatrix@generations[-1], automatonLut)
        generations.append(nextGen)
    return generations


def executeLCACommand(nCells, nGenerations, occupiedStartIndices, automatonLut):
    """
    Executes LCA command for given parameters, printing the calculated result

    arguments:
    nCells -- the size of each generation
    nGenerations -- the amount of generations to display, including starting generation
    occupiedStartIndices -- a list of indices which are occupied at start
    automatonLut -- a lookup table which describes automaton behavior

    returns:
    generations -- a list of calculated generations
    """

    # Create the first generation as a set of zero value integers
    firstGeneration = np.zeros(nCells, dtype=int)
    # Occupy the designated starting indices
    for index in occupiedStartIndices:
        firstGeneration[index - 1] = 1
    # Calculate the following generations for the specified amount minus the starting generation
    generations = propagateLCA(firstGeneration, (nGenerations-1), automatonLut)
    for generation in generations:
        printGeneration(generation)
    return generations


def parseLCACommandStringInitValues(commandString):
    """
    Gets LCA command string init values and removes initialization from the base string

    arguments:
    commandString -- a string with an LCA command

    returns:
    initValues -- a list with occupied cell indexes parsed from the command, or -1 for syntax error
    strippedString -- the command string stripped from the initializer parts, or not returned for syntax error
    """
    
    initStartIndex = commandString.find("init_start")
    initEndIndex = commandString.find("init_end")
    # Check if in the init delimiters are present, return -1 if not
    if(initStartIndex < 0):
        return -1
    if(initEndIndex < 0):
        return -1
    # Slice the initialization part from the command string
    initString = commandString[(initStartIndex + len("init_start")):initEndIndex]
    # Find the init values and map them to integers
    initValues = list(map(int, initString.split()))
    # Create a string which is stripped from initialization parts
    strippedString = commandString[:initStartIndex] + commandString[(initEndIndex + len("init_end")):]
    return initValues, strippedString


def parseLCACommandString(commandString):
    """
    Calculates and prints and LCA with behaviour and parameters specified by a command string

    arguments:
    commandString -- a command string with all parameters required for LCA display

    returns:
    generations -- a list of all printed generations based on the command

    The syntax of the command string is as follows: "aT nC nG init_start iL init_end aL"
    where:
    aT is the automaton type, corresponding to A, B or U. This argument is case sensitive.
    nC is a positive integer specifying the amount of cells in each automaton.
    nG is a positive integer specifying the amount of generations to compute and print.
    init_start is a case sensitive delimiter.
    iL is a series of whitespace-separated positive integers specifying occupied starting cell indices.
    init_end is a case sensitive delimiter.
    aL is a series of eigth whitespace-separted integer ones or zeros which define the behavior of a type U automaton.
    All command arguments are separated by one or more white spaces.
    A command with syntax errors is not guaranteed to generate a valid output.
    """

    # Parse the occupied starting indices and get the command string stripped from initialization parts
    occupiedStartIndices, strippedCommandString = parseLCACommandStringInitValues(commandString)
    # Split the command args separated by split() default whitespace
    cArgs = strippedCommandString.split()
    automatonType = cArgs[0]
    nCells = int(cArgs[1])
    nGenerations = int(cArgs[2])
    # Check the automaton type and use hard-coded values for automata type A and B or command input for type C
    if(automatonType == "A"):
        automatonLut = [0, 1, 0, 1, 1, 1, 1, 0]
    elif(automatonType == "B"):
        automatonLut = [0, 1, 1, 0, 1, 0, 1, 0]
    elif(automatonType == "U"):
        automatonLut = list(map(int, cArgs[3:]))
    # Execute the LCA with the parsed command arguments
    generations = executeLCACommand(nCells, nGenerations, occupiedStartIndices, automatonLut)
    return generations


# Run an automaton which makes a cool triangle pattern
parseLCACommandString("U 50 50 init_start 5 15 35 45 init_end 1 1 0 0 0 0 1 1")
