export AbstractBackend,CPUBackend,CPUParallelBackend,GPUBackend,OpenCLBackend,CUDABackend,MetalBackend,AMDBackend

abstract type AbstractBackend end

abstract type CPUBackend <: AbstractBackend end

abstract type CPUParallelBackend <: AbstractBackend end

abstract type GPUBackend <: AbstractBackend end

abstract type CUDABackend <: GPUBackend end

abstract type AMDBackend <: GPUBackend end