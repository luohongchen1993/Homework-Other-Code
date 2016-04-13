
# coding: utf-8

# In[1]:

f = 0.000001
sum = 0
for i in range(1000000):
    sum = sum + f
    
print sum


# In[8]:

f = 1
scale = 1000000
sum = 0
for i in range(1000000):
    sum = sum + f

intpart = sum/scale
decpart = sum%scale
decpart = decpart*10
result = "%d" %intpart
result = result+'.'

while (scale != 1):
    tmp = decpart/scale
    tmpstr = "%d" %tmp
    result = result+tmpstr
    decpart = decpart % scale
    scale = scale / 10

print result


# In[15]:

def test(numstr):
    pointflag = False
    n = len(numstr)
    s = ""
    ns = 0
    l = []
    for i in range(n):
        s=s+numstr[i]
        ns = ns+int(numstr[i])
        l.append(int(numstr[i]))
    #    print int(numstr[i]) # convert digit char to int
    #print s
    #print ns
    #print l
    


# In[28]:

def parse_str(numstr):
    if (numstr[0] == '-'):
        sign = -1
        pointflag = False
        l = []
        n = len(numstr)
        for i in range(1,n):
            if (numstr[i]=='.'):
                pointflag = True
                scale = n-i-1
                continue
            l.append(int(numstr[i]))
        if pointflag==False:
            scale = 0
        print l
        print sign
        print scale
    else:
        sign = 1
        pointflag = False
        l = []
        n = len(numstr)
        for i in range(n):
            if (numstr[i]=='.'):
                pointflag = True
                scale = n-i-1
                continue
            l.append(int(numstr[i]))
        if pointflag==False:
            scale = 0
        print l
        print sign
        print scale
    
    nstr = ""
    if (sign == -1):
        nstr += '-'
    for i in range(len(l)):
        if (i+scale==len(l)):
            nstr+='.'
        nstr+=str(l[i])
        
    print nstr
        
        
            
parse_str("-0.123456789")
parse_str("0.123456789")
parse_str("123456789")
parse_str("-12345.6789")



# In[82]:

import ipdb
class Fix:
    precision = 25
    @staticmethod
    def set_precision(p):
        Fix.precision = p
        
    def __init__(self, numstr="0.0"):
        if (numstr[0] == '-'):
            self.sign = -1
        
            pointflag = False
            self.l = []
            n = len(numstr)
            for i in range(1,n):
                if (numstr[i]=='.'):
                    pointflag = True
                    self.scale = n-i-1
                    continue
                self.l.append(int(numstr[i]))

            if pointflag==False:
                self.scale = 0
        
        else:
            self.sign = 1
            pointflag = False
            self.l = []
            n = len(numstr)
            for i in range(n):
                if (numstr[i]=='.'):
                    pointflag = True
                    self.scale = n-i-1
                    continue
                self.l.append(int(numstr[i]))
            if pointflag==False:
                self.scale = 0
    
    def __str__(self):
        nstr = ""
        zero = True
        if (self.sign == -1):
            nstr += '-'
        for i in range(len(self.l)):
            if (self.l[i] == 0 and zero == True and (i + self.scale) <= len(self.l)-2):
                continue
            if (self.l[i]!=0):
                zero = False
            if (i+self.scale==len(self.l)):
                nstr+='.'
            nstr+=str(self.l[i])
        return nstr
    
    
    def printnum(self):
        nstr = self.__str__()
        print nstr
            
    def change_scale(self, newscale):
        it = self.scale
        while (it<newscale):
            self.l.append(0)
            it = it+1
        self.scale = newscale
    
    #def simplify(self):
    #    print "erase some zeros in front"
    
    def cmp_abs_val(self,other):
        max_scale = max(self.scale,other.scale)
        self.change_scale(max_scale)
        other.change_scale(max_scale)
        i1 = 0
        i2 = 0
        n1 = len(self.l)
        n2 = len(other.l)
        while(i1<n1 and self.l[i1]==0):
            i1 = i1+1
            if i1 == n1: break
        
        while(i2<n2 and other.l[i2]==0):
            i2 = i2+1
            if i2 == n2: break
                
        if (n1-i1>n2-i2):
            return 1
        elif n1-i1<n2-i2:
            return -1
        else:
            while (i1<n1 and i2<n2):
                if self.l[i1]>other.l[i2]:
                    return 1
                elif self.l[i1]<other.l[i2]:
                    return -1
                else:
                    i1=i1+1
                    i2=i2+1
            return 0
        
    def __gt__(self,other):
        if (self.sign==1 and other.sign==-1):
            return True
        elif (self.sign==-1 and other.sign==1):
            return False
        elif (self.sign==1 and other.sign==1):
            if (self.cmp_abs_val(other)==1):
                return True
            else:
                return False
        else:
            if (self.cmp_abs_val(other)==-1):
                return True
            else:
                return False
    
    def __lt__(self,other):
        if (self.sign==1 and other.sign==-1):
            return False
        elif (self.sign==-1 and other.sign==1):
            return True
        elif (self.sign==1 and other.sign==1):
            if (self.cmp_abs_val(other)==-1):
                return True
            else:
                return False
        else:
            if (self.cmp_abs_val(other)==1):
                return True
            else:
                return False
    
    def __eq__(self, other):
        return ((not self>other) and (not self<other))
    
    def __ne__(self, other):
        return not self == other

    
    def __ge__(self, other):
        return (self > other) or (self == other)
    
    def __le__(self, other):
        return (self < other) or (self == other)
    
    def __neg__(self):
        newins = self
        newins.sign = (-1)*self.sign
        return newins
    
    def __abs__(self):
        newins = self
        newins.sign = 1
        return newins
    
    def __add__(self,other):
        result = Fix("0.0")
        
        if (self.sign==other.sign):
            max_scale = max(self.scale,other.scale)
            self.change_scale(max_scale)
            other.change_scale(max_scale)
            n1 = len(self.l)
            n2 = len(other.l)
            
            i1 = n1-1
            i2 = n2-1
            numvec = []
            carry = 0
            while (i1>=0 and i2>=0):
                tmp = self.l[i1]+other.l[i2]+carry
                res = tmp % 10
                carry = tmp / 10
                numvec.insert(0,res)
                i1=i1-1
                i2=i2-1
            
            while (i1>=0):
                tmp = self.l[i1]+carry
                res = tmp % 10
                carry = tmp / 10
                numvec.insert(0,res)
                i1 = i1-1
            
            while (i2>=0):
                tmp = other.l[i2]+carry
                res = tmp % 10
                carry = tmp / 10
                numvec.insert(0,res)
                i2=i2-1
            
            result.l = numvec
            result.scale = max_scale
            result.sign = self.sign
            return result
        else:
            max_scale = max(self.scale,other.scale)
            self.change_scale(max_scale)
            other.change_scale(max_scale)
            if self.cmp_abs_val(other)==1:
                num1 = Fix(str(self))
                num2 = Fix(str(other))
            elif self.cmp_abs_val(other)==-1:
                num1 = Fix(str(other))
                num2 = Fix(str(self))
            else:
                return Fix("0.0")
            
            n1 = len(num1.l)
            n2 = len(num2.l)
            
            i1 = n1-1
            i2 = n2-1
            numvec = []
            carry = 0
            while (i1>=0 and i2>=0):
                if (num1.l[i1]+carry>=num2.l[i2]):
                    tmp = num1.l[i1]+carry-num2.l[i2]
                    carry = 0
                    numvec.insert(0,tmp)
                else:
                    tmp = num1.l[i1]+carry-num2.l[i2]+10
                    carry = -1
                    numvec.insert(0,tmp)

                i1=i1-1
                i2=i2-1
            
            while (i1>=0):
                if (num1.l[i1]+carry>=0):
                    tmp = num1.l[i1]+carry
                    carry = 0
                    numvec.insert(0,tmp)
                else:
                    tmp = num1.l[i1]+carry+10
                    carry = -1
                    numvec.insert(0,tmp)

                i1 = i1-1
            
            result.l = numvec
            result.scale = max_scale
            neg_other = -num2
            if(num1>=neg_other):
                result.sign = 1
            else:
                result.sign = -1
            
            return result
    
    def __sub__(self,other):
        neg_other = Fix(str(other))
        neg_other.sign = (-1)*neg_other.sign
        result = self + neg_other
        return result
    
    def times10(self):
        self.l.append(0)
        
    def divide10(self):
        self.l.insert(0,0)
        self.scale = self.scale + 1
    
    def __mul__(self,other):
        tmp = Fix(str(self))
        num = other.l
        memory = Fix("0.0")
        for i in range(len(num)):
            curr_num = num[i]
            
            if curr_num==0:
                if (i != len(num) - 1):
                    memory.times10()
                continue
            else:
                for j in range(0,curr_num):
                    memory = memory + tmp
                                        
                if (i != len(num) - 1):
                    memory.times10()
                    
        for i in range(other.scale):
            memory.l.insert(0,0)
        
        memory.scale = memory.scale+other.scale
        memory.sign = self.sign*other.sign
        return memory
    
    def __div__(self,other):
        dividend = Fix(str(self))
        num = self.l
        digit_count = 0
        num_vec = []
        multiplier = Fix("1.0")
        dividend = abs(dividend)
        f = abs(other)
        while (dividend>f): # end with f>=div
            f.times10()
            multiplier.times10()
        
        if (dividend==f): # f==div
            result = multiplier
            result.sign = self.sign*other.sign
            return result
        
        # f>div
        while (f>dividend):
            dividend.times10();
            multiplier.divide10();
            digit_count+=1
            if (digit_count == Fix.precision):
                break

        # either digit count ==  precision or div>=f
        # end with div>f but div<10*f
        
        div = dividend
        result_scale = -1
        zero_const = Fix("0.0")
        while (digit_count<Fix.precision):
                curr_num = 0;
                
                while (div>=f):
                    div = div-f
                    curr_num = curr_num+1
                num_vec.append(curr_num)
                result_scale+=1
                digit_count+=1
                if (div==zero_const):
                    break
                div.times10()

                    
        if (result_scale == -1):
            return Fix("0.0")
        else:
            result = Fix("0.0")
            result.l = num_vec
            result.scale = result_scale
            result.sign = self.sign*other.sign
            result = result*multiplier
            return result


# In[87]:

a = Fix("-12.345")
b = Fix("7.98")
a = a/b
a.printnum()
#result = a/b
#result.printnum()
#result.printnum()
#a.printnum()
#b.printnum()
#c=-b
#c.printnum()
#print a.l
#print a.sign
#print b.sign
#print (a<=b)
#print str(a)

#a.change_scale(10)
#a.printnum()
#print b.scale
#print a.cmp_abs_val(b)
#print a.to_string()

