
struct part_struct
     coorx::Float64
     coory::Float64
     coorz::Float64
     velx::Float64
     vely::Float64
     velz::Float64
     gradx::Vector{Float64}
     grady::Vector{Float64}
     gradz::Vector{Float64}
end 

#create empty array
num_res = 10
particle = Array{part_struct,1}(undef, num_res)

##initialize array
for i = 1:num_res
      particle[i] = part_struct( rand(), rand(), rand(), rand(), rand(), rand(), [rand(); rand()], [rand(); rand()], [rand(); rand()] )
end

for i = 1:num_res
     println("information ",particle[i] )
end
#
#println("which ", particle[2].coorx, particle[2].gradx[1], particle[2].gradx[2])
