import math

def function(x):
    #return (x-math.sqrt(5))*(x+math.sqrt(5))
    #return (math.sin(x+1))
    return math.log(x-math.sqrt(15))


if __name__ == "__main__":
    x = 0
    a = 4 #left x value
    b = 20 # right x value
    last_middle_value = ""

    leftAbove = function(a) > 0
    rightAbove = function(b) > 0
    go_on = True

    # We check if either bound is already a zero:
    if (function(a) == 0):
         go_on = False
         middle_value = a

    if (function(b) == 0):
        go_on = False
        middle_value = b

    # we check if the specified a, b and function satisfy the condition.
    if (go_on and rightAbove and leftAbove or not rightAbove and not leftAbove):
        go_on = False
        print("invalid a,b or function specified")
    


    while go_on:
        middle_value = (b+a)/2
        print("current value:", middle_value)
        result = function(middle_value)
        print("range: from", a, "to", b)
        
        if leftAbove:
            if (result > 0):
                a = middle_value
            elif (result < 0):
                b = middle_value
        else: # rightabove
            if (result < 0):
                a = middle_value
            elif (result > 0):
                b = middle_value
       
        if (last_middle_value == middle_value): 
            break

        last_middle_value = middle_value
    
    print("one zero is:", middle_value)