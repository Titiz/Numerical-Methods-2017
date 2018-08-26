

def add(string1, string2, length):
    i = length - 1
    result = ""
    remainder = 0
    while (i> -1):
        term = 0
        a = string1[i]
        b = string2[i]
        if remainder == 1:
            term += 1
            remainder = 0
        if a == "1":
            term += 1
        if b == "1":
            term += 1
        if term == 0:
            result += "0"
        elif term == 1:
            result += "1"
        elif term == 2:
            result += "0"
            remainder += 1
        elif term == 3:
            result += "1"
            remainder += 1
        i -= 1
        print("a:", a, "b:", b)
        print("result:", result[::-1])
        print("remainder:", remainder)
    if (remainder == 1):
        result += "1"
    print("length of result", len(result), ":", length)
    return result[::-1]


s1 = "11001"
s2 = "01010"

string1 = "0111111101111110110001010011101010001101010010010001"
string2 = "0111111101111110110001010011101010001101010010010101"
            
        
print(add(string1, string2, len(string1)))
print(add(s1, s2, len(s1)))

