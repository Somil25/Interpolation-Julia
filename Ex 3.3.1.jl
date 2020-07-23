
# Excercise 3.3.1
# Polynomial Interpolation (Newton)
function polyinter(x::Array,y::Array,z)
         n=length(x) # n is the number of data points & degree of the polynomial is n-1
         a=Array{Float64}(undef,n)
         for i in 1:n
             a[i]=y[i]
         end
         for j in 2:n
            for i in reverse(collect(j:n)) # reverse function is used to reverse order from n to j 
                a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)]) # keeps on updating the differences
            end 
        end
        l=a # divided difference
        sum=l[1]
        prefac=1.0
        for k in 1:(n-1)
             prefac=prefac*(z-x[k])
             sum=sum+l[k+1]*prefac # f[x0] + f[x0; x1](z - x0) + ... + f[x0; x1 ;...; xn](z - x0)(z -x1)...(z - xn-1) (Newton)
        end
        sum
end


# Hermite Interpolation Polynomial
function hermpoly(x::Array,y::Array,yp::Array,p)
         n=length(x) # n number of data points.
         l=2n
         z=Array{Float64}(undef,l)
         a=Array{Float64}(undef,l)
         for i in 1:n # Since z2i = z2i-1 = xi
             z[2i-1]=x[i]
             z[2i]=x[i]
         end
         for i in 1:n
             a[2i-1]=y[i]
             a[2i]=y[i]
         end
         for i in reverse(collect(2:n)) # computing the first divide differences 
             a[2i]=yp[i]
             a[2i-1]=(a[2i-1]-a[2i-2])/(z[2i-1]-z[2i-2])
         end
         a[2]=yp[1]
         for j in 3:l #computing higher order divided differences starting from 3 since a(1) & a(2) are already updated.
             for i in reverse(collect(j:l))
                 a[i]=(a[i]-a[i-1])/(z[i]-z[i-(j-1)])
              end
         end
         sum=a[1]
         prefac=1.0
         for k in 1:(2n-1) #H_n(x) = f[z0] +Σ_[i=1]^[n] f[z0,.., zi](x -z_(0)) ....(x -z_(i-1)) (Hermit)
             prefac=prefac*(p-z[k])
             sum=sum+a[k+1]*prefac 
         end
         return sum
end

#Plotting the two interpolating polynomials together with f(x)
using PyPlot  # make sure to use before plotting
f(x)=exp(x)+sin(10*x) # function defined
xi=0:1/100:3 
#data
x=[0,0.4,1 ,2, 2.6 ,3]
y=[1,0.735,2.17,8.30,14.2,19.1]
yp=[11,-5.04,-5.67,11.5, 9.9,21.6]
funct=map(f,xi)
d=map(w->hermpoly(x,y,yp,w),xi)
pol=map(z->polyinter(x,y,z),xi)
plot(xi,pol,label="Polynomial interpolant")
plot(xi,d,label="Hermite interpolation")
plot(xi, funct, label="f(x)")
scatter(x, y, label="data")
xlabel("X-axis")
ylabel("Y-axis")
legend(loc="upper right");
