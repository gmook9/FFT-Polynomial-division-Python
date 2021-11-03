#Katie Pell and Garret Mook
#Dr. Ravikumar - CS415
#Fall 2021

import math
from math import ceil, floor

#Get elements at even indices
def getEvens(A):
    if len(A) <= 2:
        return [A[0]]
    E = []
    for i in range(len(A)):
        if i % 2 == 0:
            E.append(A[i])
    return E

#Get elements at odd indices
def getOdds(A):
    if len(A) <= 2:
        return [A[1]]
    O = []
    for i in range(len(A)):
        if i % 2 != 0:
            O.append(A[i])

    return O

#Function to do rounding. Different from clean which removes extra 0s from end.
def rounding(A):
    for i in range(len(A)):
        x = A[i].real
        y = A[i].imag
        newx = x
        newy = y
        if abs(floor(x) - x) <= 0.000001 and abs(floor(x) - x) != 0:
            newx = floor(x)
            newy = y
        elif abs(ceil(x) - x) <= 0.000001 and abs(ceil(x) - x) != 0:
            newx = ceil(x)
        if abs(floor(y) - y) <= 0.000001 and abs(floor(y) - y) != 0:
            newy = floor(y)
        elif abs(ceil(y) - y) <= 0.000001 and abs(ceil(y) - y) != 0:
            newy = ceil(y)

        if x < 0 and int(x) - x <= 0.000001 or int(x) + 1 - x <= 0.000001:
            newx = int(newx)

        if newy == 0:
            A[i] = round(newx, 6)
        else:
            A[i] = complex(round(newx, 6), round(newy, 6))

    return A

def padBoth(u, v):
    while len(u) < len(v):
        u.append(0)
    while len(v) < len(u):
        v.append(0)

def pad(v):
    x = []
    for item in v:
        x.append(item)
    m = 1
    #Make sure the list is nearest power of 2 equal or greater.
    while m < len(x):
        x.append(0)
        m = m * 2
    return x

#Simple function to do component wise multiplication
#*****************EXAMPLE*****************
# u = [1,2,3,4]
# v = [5,6,7,8]
# c = [(1*5),(2*6),(3*7),(4*8)]
# Therefore:
# c = [5,12,21,32]
#*****************************************
def componentMultiply(u, v):
    c = [0] * len(u)
    for i in range(len(u)):
        c[i] = (u[i]*v[i])
    
    return c

#Takes in a list as a parameter which will be plugged into FFT
#function along with the parameter W that is computer here. This exists
#as a way to not have to enter W each time we want to do a FFT. We just need
#one parameter and thats the list.
def FFT1(A):
    #r is going to be the W parameter that we pass into FFT. 
    r = FFT(A,complex(math.cos(2*math.pi/len(A)), math.sin(2*math.pi/len(A))))
    #Save the first element since it will be popped
    temp = r[0]
    #Remove the element at index 0
    r.pop(0) 
    #Reverse the list
    r.reverse()
    #Reinsert the first element back into the first position
    r.insert(0, temp)
    return r

def FFT(A, w):
    if len(A)==1:
        return A
    Ae = FFT(getEvens(A), w*w)
    Ao = FFT(getOdds(A), w*w)
    r = []
    for k in range(len(Ae)):
        # if k < len(getEvens(A)) and k < len(getOdds(A)):
        r.append(Ae[k] + (w ** k) * Ao[k])
    for k in range(len(Ae)):
        # if k < len(getEvens(A)) and k < len(getOdds(A)):
        r.append(Ae[k] - (w ** k) * Ao[k])
    return r

def IFFT(A):
    x = FFT1(A)
    for j in range(len(A)): # for every element
        x[j] = x[j]/len(A) # return that element, normalized
    #Save the first element since it will be popped
    temp = x[0]
    #Remove the first element
    x.pop(0)
    #Reverse the list
    x.reverse() 
    #Reinsert the first element back into the first position
    x.insert(0, temp) 

    return clean(x)

def Convolution(u, v):
    #Initialize a list to append the items from u and v into
    w = []
    x = []
    for item in u:
        w.append(item)
    for item in v:
        x.append(item)
    curr = 1
    #Make sure the list is nearest power of 2 equal or greater. This is the same
    #as what we do in the pad function
    while curr < len(x):
        x.append(0)
        curr = curr * 2
    #pads w and x
    while len(w) < curr:
        w.append(0)
    while len(x) < curr:
        x.append(0)

    # now w and x are padded correctly
    w = FFT1(w)
    x = FFT1(x)

    #do component wise multiplication on the elements
    y = componentMultiply(w, x)
    #Then take the IFFT of the result and be sure to apply rounding
    y = rounding(IFFT(y))

    #This works like clean, pops the last element if its equal to 0
    while y[-1] == 0:
        y.pop(-1)

    return y

#Function to clean. If the last element of v is 0
#we remove it with pop. pop takes in a index. We give it
#the index of -1 because that means the last element.
#*****************EXAMPLE*****************
# v = [1,2,3,4,0,0,0]
# Since its a while loop its pops everytime the last
# element is still equal to 0. After the first pop we get:
# v = [1,2,3,4,0,0,0]
# The while statement is still true so we pop again the last element
# v = [1,2,3,4,0,0]
# we repeat this until the last element is not zero
# v = [1,2,3,4,0]
# The v that is returned then is:
# v = [1,2,3,4]
#****************************************
def clean(v):
    while v[-1] == 0:
        v.pop(-1)
    return v

#Function to find the coefficients of a polynomial given the roots.
def Poly(v):
    #Base case
    if len(v) == 1:
        v[0] = -v[0]
        v.append(1)
        return v
    w = v
    #Python is not inclusive with ranges for the right side of :, so we must go up to len(w)//2
    #However, it the element at len(w)//2 is NOT included in the list v1
    #instead that element is picked up by v2 when it is on the left side of :
    #*****************EXAMPLE*****************
    #w = [1,2,3,4,5,6]
    #v1 = w[0 : len(w)//2], v1 = [1,2,3]
    #v2 = w[len(w)//2: len(w)], v2 = [4,5,6]
    #****************************************
    v1 = w[0 : len(w)//2]
    v2 = w[len(w)//2: len(w)]
    #Recursive Calls
    u1 = Poly(v1)
    u2 = Poly(v2)
    u = Convolution(u1, u2)
    u = clean(u)
    return u

#Helper function to find the difference from the elements at each index
#between each of the two lists.
#*****************EXAMPLE*****************
# a = [5,6,7,8]
# b = [1,2,3,4]
# c = divideHelper(combineBothListsForLoop)
# c = [(5-1),(6-2),(7-3),(8-4)]
# c = [(4),(4),(4),(4)]
#*****************************************
def divideHelper(combineBothListsForLoop):
    firstFile = [firstFileElementIdx - secondFileElementIdx for firstFileElementIdx, secondFileElementIdx in combineBothListsForLoop]
    return firstFile

#Helper function for polynomial Division
#Function used to combine lists.
#*****************EXAMPLE*****************
# a = [1,2,3,4]
# b = [5,6,7,8]
# c = combineLists(a, b)
# c = [(1,5),(2,6),(3,7),(4,8)]
#*****************************************
def combineLists(firstList, quotientFoundTimesSecondList): 
	totalListSize = len(firstList) if len(firstList) < len(quotientFoundTimesSecondList) else len(quotientFoundTimesSecondList) 
	combinedList = [] 
	for i in range(totalListSize): 
		combinedList.append((firstList[i], quotientFoundTimesSecondList[i])) 
	return combinedList

#Helper for divide function
def divideHelperB(quotientFound,secondFile):
    quotientFoundTimesSecondList = [quotientFound * elementIdx for elementIdx in secondFile]
    return quotientFoundTimesSecondList

# This is the pseudocode for the algorithm that we used for divide (different from the one given to us)
    # function polyDivide(num, den):
    # input: lists of size m and n respectively where they represent
    #     polynomials A(x) and B(x) with those coefficients
    # output: list with coefficients of the polynomial C(x) where C(x) = A(x) / B(x)

    # Remove trailing zeroes from num and den

    # If it's empty, add a zero

    # If the numerator has degree larger than the denominator, add zeroes until they are the same

    # Create list that holds result; r = []

    # Divisor = Coefficient of the highest-degree term in the denominator, i.e. den(n - 1)

    # For i in range (number of zeroes added to denominator + 1):
    #     Calculate ratio = numerator first coefficient / denominator first coefficient

    #     Add ratio to list of results

    #     If ratio is not zero:
    #         Create a copy of denominator and multiply it by ratio; d = []
    #         Create a list e such that:
    #             e = {{num[0],d[0]},{num[1],d[1]},...,{num[n-1],d[n-1]}}
    #         Add to the result list {a-b} for {a,b} in e

    #         Remove the first element from the numerator
    #         Remove the second element from the denominator

    # Trim trailing zeroes
def Divide(firstFile, secondFile):
    #Initialize a list for the quotient 
    quotientList = []
    #Get the lengths of the lists from each file
    lengthOfFirstFileList = len(firstFile)
    lengthOfSecondFileList = len(secondFile)
    differenceInListSize = lengthOfFirstFileList - lengthOfSecondFileList
    #Pad if needed
    #*****************EXAMPLE*****************
    #If firstFile = [1 2 3 4 5]
    #and firstFile =[6 7]
    #Then after shifting to match degrees the 
    #firstFile =  [1 2 3 4 5]
    #secondFile = [0 0 0 6 7]
    #****************************************
    padWithZero = [0]
    secondFilePadded = padWithZero * differenceInListSize
    secondFile = secondFilePadded + secondFile
    #if we look at the list from our second file
    #we want to now grab the LAST element and make it a factor
    #or instead of factor your could think of this as what we
    #are going to divide by
    #*****************EXAMPLE*****************
    #secondFile = [0 0 0 6 7]
    #therefore, factor would be 7
    #factor is equivalent to secondFile[len(secondFile)-1]
    #****************************************
    #No k++ needed in the loop because for loop iterates automatically
    for k in range(differenceInListSize + 1):
        #Divide the LAST element from the first list by the last 
        #element of the second list.
        lastIdxOfFirstFile = len(firstFile)-1
        lastIdxOfSecondFile = len(secondFile)-1
        quotientFound = firstFile[lastIdxOfFirstFile] / secondFile[lastIdxOfSecondFile]
        #This becomes the quotient. We create a list of quotients so each time we find
        #a quotient it goes into the list of quotients. This continues until
        #quotient is a list of all the quotients found. This is the list that will
        #ultimately hold our answer
        quotientList = [quotientFound] + quotientList
        zeroConstant = 0
        #When this if runs it will call a helper which will 
        #aid by subtracting quotientFound * secondFile from firstFile 
        #unless the quotientFound is equal to zero 0
        if quotientFound != zeroConstant:
            #quotient found * 2ndList at x (loops through all elements in list)
            quotientFoundTimesSecondList  = divideHelperB(quotientFound,secondFile)
            #then take that answers and subtract it from elements in the firstlist
            combineBothListsForLoop = combineLists(firstFile, quotientFoundTimesSecondList)
            firstFile = divideHelper(combineBothListsForLoop)
            # firstFile = [firstFileElementIdx - secondFileElementIdx for firstFileElementIdx, secondFileElementIdx in combineBothListsForLoop]
        #remove LAST element of firstList
        firstFile.pop()
        #remove first element of 2ndList
        secondFile.pop(zeroConstant)
    #once the for loops completes we have compiled a full list of all the quotients
    #this list is the answer.
    return quotientList

#-----------------Divide from Dr.Ravikumar Pseudocode------------------------
#We spoke in office hours about using our new divide algorithm above and got the okay. Here is our old one that we formatted based of your pseudocode
#that way you can see we did both. Thank you.
#
# def Divide(a, b):
#     # a is a vector of length n
#     # b is a vector of length m
#     # m <= n
#     # a and b are polynomials where a(x) = b(x) * c(x)
#     # returns c(x) as a vector
#     a = pad(a)
#     b = pad(b)
#     padBoth(a, b)
#     f1 = FFT1(a)
#     f2 = FFT1(b)
#     f3 = [] # creates a vector of same size as f1 but all zeroes
#     for item in f1:
#         f3.append(item)
#     for item in f1:
#         item = 0
#     for j in range(len(f2)): # check every element in f2
#         if f2[j] == 0: # correct, as required by problem
#             # for item in f1: 
#             #     item += 0.0000001 * len(f1) # add a small coefficient to each x in b
#             f3[j] = f1[j] / (f2[j] + 0.00000001)
#         else:
#             f3[j] = f1[j] / f2[j]
        
#     f4 = IFFT(f3)
#     return rounding(f4)
#
#----------------------------------------------------------------------------

#Function that will be called to print the menu
def printMenu():
    print("1. Problem 1")
    print("2. Problem 2")
    print("3. Quit")

#Function that will be called to create a list [] from
#a specific file given a file name passed in as a parameter
def vectorFromFile(filename):
    fileA = open(filename, "r")
    fileContents = fileA.read()
    fileVector = fileContents.split(" ")
    for i in range(len(fileVector)):
        fileVector[i] = float(fileVector[i])
    fileA.close()
    return fileVector

#Function that will control the user menu
def userMenu():
    printMenu()
    answer = int(input("Enter an option from the menu above: "))
    while answer != 3:
        if answer == 1:
            filename1 = input("Enter filename: ")
            v = vectorFromFile(filename1) # Open and read the file, return vector with each item in it
            print(Poly(v)) # Problem 1 answer
        elif answer == 2:
            filename2 = input("Enter first filename: ") # Get first vector
            w = vectorFromFile(filename2)
            p = w[:]
            filename3 = input("Enter second filename: ") # Get second vector
            x = vectorFromFile(filename3)
            q = x[:]
            print(Divide(p, q)) # Problem 2 answer
        else:
            answer = int(input("Enter a valid option from the menu above: "))
            continue
        printMenu()
        answer = int(input("Enter an option from the menu above: "))

def main():
    userMenu()
   
main()    