# Module.jl is already loaded by Preprocessing_routines.jl (included from Main.jl).
# Just bring G into scope here.
using .CFDVars: G
using Printf

function run_output_dir()
    if G.exec_mode == 0
        return "serial"
    elseif G.exec_mode == 1
        if G.parallel_mode == 1
            return "parallel_1"
        elseif G.parallel_mode == 2
            return "parallel_2"
        elseif G.parallel_mode == 3
            return "parallel_3"
        elseif G.parallel_mode == 4
            return "parallel_4"
        else
            return "parallel_unknown"
        end
    end
    return "run_unknown"
end

# Fortran unformatted binary record writer
function write_record(io, data)
    n = sizeof(data)
    write(io, Int32(n))
    write(io, data)
    write(io, Int32(n))
end

function write_record(io, val::Int32)
    write(io, Int32(4))
    write(io, val)
    write(io, Int32(4))
end

# Write grid, flow, and restart files (binary, Fortran-compatible).
function output!(flag)

    nvars = G.nprims + G.nconserv + 9

    # grid file
    mkpath(run_output_dir())
    open(joinpath(run_output_dir(), "grid.xyz"), "w") do fgrid
        write_record(fgrid, Int32(G.nblocks))
        dims = Int32[]
        for nbl = 1:G.nblocks
            push!(dims, Int32(G.NI[nbl]))
            push!(dims, Int32(G.NJ[nbl]))
            push!(dims, Int32(G.NK[nbl]))
        end
        write_record(fgrid, dims)
        for nbl = 1:G.nblocks
            coords = vcat(vec(G.xgrid[:,:,:,nbl]),
                          vec(G.ygrid[:,:,:,nbl]),
                          vec(G.zgrid[:,:,:,nbl]))
            write_record(fgrid, Float64.(coords))
        end
    end

    # flow file
    filename = (flag == 1) ? @sprintf("flow%05d.xyz", G.iter) : "flow.xyz"
    flow_path = joinpath(run_output_dir(), filename)

    open(flow_path, "a") do fflow
        write_record(fflow, Int32(G.nblocks))
        dims = Int32[]
        for nbl = 1:G.nblocks
            push!(dims, Int32(G.NI[nbl]))
            push!(dims, Int32(G.NJ[nbl]))
            push!(dims, Int32(G.NK[nbl]))
            push!(dims, Int32(nvars))
        end
        write_record(fflow, dims)
        for nbl = 1:G.nblocks
            data = Float64[]
            for n = 1:G.nprims;   append!(data, vec(G.Qp[:,:,:,nbl,n])) end
            for n = 1:G.nconserv; append!(data, vec(G.Qc[:,:,:,nbl,n])) end
            append!(data, vec(G.ix[:,:,:,nbl]))
            append!(data, vec(G.iy[:,:,:,nbl]))
            append!(data, vec(G.iz[:,:,:,nbl]))
            append!(data, vec(G.jx[:,:,:,nbl]))
            append!(data, vec(G.jy[:,:,:,nbl]))
            append!(data, vec(G.jz[:,:,:,nbl]))
            append!(data, vec(G.kx[:,:,:,nbl]))
            append!(data, vec(G.ky[:,:,:,nbl]))
            append!(data, vec(G.kz[:,:,:,nbl]))
            write_record(fflow, data)
        end
    end

    # restart file
    open(joinpath(run_output_dir(), "restart.xyz"), "w") do fflow
        write_record(fflow, Int32(G.nblocks))
        dims = Int32[]
        for nbl = 1:G.nblocks
            push!(dims, Int32(G.NI[nbl]))
            push!(dims, Int32(G.NJ[nbl]))
            push!(dims, Int32(G.NK[nbl]))
            push!(dims, Int32(nvars))
        end
        write_record(fflow, dims)
        for nbl = 1:G.nblocks
            data = Float64[]
            for n = 1:G.nprims;   append!(data, vec(G.Qp[:,:,:,nbl,n])) end
            for n = 1:G.nconserv; append!(data, vec(G.Qc[:,:,:,nbl,n])) end
            write_record(fflow, data)
        end
    end

    println("Output written: ", flow_path)
end

# Memory is reclaimed automatically by Julia's GC.
function deallocate_routine!()
    println("Deallocate Done (GC handles memory in Julia)")
end

function volume_integral!()
end