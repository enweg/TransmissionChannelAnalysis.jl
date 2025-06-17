@testset "_find_nonmissing_period" begin
    X = convert(Matrix{Union{Missing,Float64}}, randn(10, 3))
    idx_first, idx_last = TransmissionChannelAnalysis._find_nonmissing_period(X)
    @test idx_first == 1
    @test idx_last == 10

    X[1, 1] = missing
    idx_first, idx_last = TransmissionChannelAnalysis._find_nonmissing_period(X)
    @test idx_first == 2
    @test idx_last == 10

    X[2, 3] = missing
    X[9, 1] = missing
    idx_first, idx_last = TransmissionChannelAnalysis._find_nonmissing_period(X)
    @test idx_first == 3
    @test idx_last == 8
end

@testset "_find_data_overlap" begin
    X = convert(Matrix{Union{Missing,Float64}}, randn(10, 3))
    Y = convert(Matrix{Union{Missing,Float64}}, randn(10, 3))
    Z = convert(Matrix{Union{Missing,Float64}}, randn(10, 3))

    Xout, = TransmissionChannelAnalysis._find_data_overlap(X)
    @test Xout == X

    X[1, 1] = missing
    Xout, Yout = TransmissionChannelAnalysis._find_data_overlap(X, Y)
    @test Xout == X[2:end, :]
    @test Yout == Y[2:end, :]

    Z[1:3, :] .= missing
    Xout, Yout, Zout = TransmissionChannelAnalysis._find_data_overlap(X, Y, Z)
    @test Xout == X[4:end, :]
    @test Yout == Y[4:end, :]
    @test Zout == Z[4:end, :]

    Y[7:end, :] .= missing
    Xout, Yout, Zout = TransmissionChannelAnalysis._find_data_overlap(X, Y, Z)
    @test Xout == X[4:6, :]
    @test Yout == Y[4:6, :]
    @test Zout == Z[4:6, :]
end
