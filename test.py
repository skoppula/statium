def binary_search(elements, e):
    left, right = 0, len(elements)-1
    mid = (left+right)/2
    
    if(elements[left] == e):
        return left
    elif(elements[right] == e):
        return right

    print(left, right)

    while(right - left > 1):
	print(left, right)
        if(e == elements[mid]):
            return mid
        elif(e > elements[mid]):
            left = mid
        else:
            right = mid
            
        mid = (left+right)/2
    
    print('Could not find element in sorted list')   
    return -1

print(binary_search([-10,2,3, 3.5,4,4.5,5,6,7.5,8,9,14], 6))
