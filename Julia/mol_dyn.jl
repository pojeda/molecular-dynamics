
struct part_struct
     x::Float64
     y::Float64
     z::Float64
end 

#outer constructor
part_struct() = part_struct(rand(),rand(),rand())

#create empty array
num_res = 10
particle = Array{part_struct,1}(undef, num_res)

#initialize array
for i = 1:num_res
     particle[i] = part_struct()
end

for i = 1:num_res
     println("information",particle[i] )
end

println("second", particle[2].x, particle[2].y, particle[3].z)
