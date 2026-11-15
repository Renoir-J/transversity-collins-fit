@testset "synthetic fixed times: ownership, labels, provenance and setting rejection" begin
    # Explicit synthetic provenance; it cannot admit physical or production work.
    qs=[2.4,10.,1000.];ts=[0.,.06,.2];sha=repeat("0",64)
    profile=PreparedEvolutionTimes(qs,ts,2.4,5,.327,"synthetic","unit-test",sha)
    qs[2]=11.;ts[2]=.07
    @test profile.Q2==(2.4,10.,1000.) && profile.tau==(0.,.06,.2)
    @test !ismutabletype(typeof(profile))
    @test prepared_times(profile,[1000.,2.4,10.,10.],2.4)==[.2,0.,.06,.06]
    @test evolution_time(profile,10.,2.4)==.06
    @test prepared_time_metadata(profile).source_sha256==sha
    clone=PreparedEvolutionTimes(profile.Q2,profile.tau,2.4,5,.327,"synthetic","other-source",sha)
    @test CollinsFit.numerical_snapshot(profile)!=CollinsFit.numerical_snapshot(clone)
    clone2=PreparedEvolutionTimes(profile.Q2,profile.tau,2.4,5,.327,"other-model","unit-test",sha)
    @test CollinsFit.numerical_snapshot(profile)!=CollinsFit.numerical_snapshot(clone2)
    @test_throws Exception setfield!(profile,:nf,4)
    @test_throws Exception setindex!(profile.tau,.1,2)
    for q in (1.,nextfloat(10.),Inf,NaN,true)
        @test_throws ErrorException prepared_times(profile,[q],2.4)
    end
    @test_throws ErrorException prepared_times(profile,Float64[],2.4)
    @test_throws ErrorException prepared_times(profile,[10.],3.)
    @test_throws ErrorException prepared_times(profile,[10.],2.4;nf=4)
    @test_throws ErrorException prepared_times(profile,[10.],2.4;nf=true)
    @test_throws MethodError alpha_s(profile,10.)
    @test_throws ErrorException ForwardDiff.derivative(q->only(prepared_times(profile,[q],2.4)),10.)
    @test_throws ErrorException ForwardDiff.derivative(t->PreparedEvolutionTimes([2.4,10.],[0.,t],2.4,5,.327,"synthetic","test",sha).alpha0,.06)
    good=([2.4,10.,1000.],[0.,.06,.2],2.4,5,.327,"synthetic","unit-test",sha)
    replacements=[(1,[2.4,2.4,1000.]),(1,[10.,2.4,1000.]),(1,[0.,10.,1000.]),
        (1,[2.4,NaN,1000.]),(1,Float64[]),(2,[0.,.06]),(2,[0.,-.06,.2]),
        (2,[0.,.06,.06]),(2,[0.,.06,1.01]),(2,[.001,.06,.2]),(2,[0.,NaN,.2]),
        (3,3.),(3,true),(4,true),(4,7),(4,5.),(5,0.),(5,NaN),(5,1.1),
        (6,""),(6,"space invalid"),(7,""),(7,repeat("a",257)),(8,"a"),(8,repeat("A",64))]
    for (index,value) in replacements
        args=Any[good...];args[index]=value
        @test_throws ErrorException PreparedEvolutionTimes(args...)
    end
    backward=PreparedEvolutionTimes([1.,2.4,10.],[-.03,0.,.06],2.4,5,.327,"synthetic","unit-test",sha)
    @test evolution_time(backward,1.,2.4)==-.03
    @test_throws ErrorException HoppetPointPlan("absent",[.1],[10.],2.4,profile)
    @test_throws ErrorException HoppetMatchingPlan("absent",[.1],[10.],2.4,profile;process=:SIDIS)
end

@testset "synthetic prepared dense operator: default identity, AD and mutation guards" begin
    coupling=AsymptoticTwoLoopCoupling();qs=[2.4,10.,1000.]
    times=evolution_time.(Ref(coupling),qs,2.4)
    profile=PreparedEvolutionTimes(qs,times,2.4,5,alpha_s(coupling,2.4),"synthetic-asymptotic","unit-test",repeat("0",64))
    grid=collect(range(.001,1.;length=31));h=grid.*(1 .-grid).^3
    old=collinear_lo_operator(grid,10.,2.4,coupling)
    op=collinear_lo_operator(grid,10.,2.4,profile)
    @test old.time_profile===nothing && old.prepared_signature===nothing
    @test old.coupling_model==:asymptotic_two_loop && old.lambda_GeV==.225
    @test old.matrix==op.matrix
    @test apply_operator(old,h)==apply_operator(op,h)
    @test op.time_profile===profile && op.coupling_model==:prepared_times && op.lambda_GeV===nothing
    grid[2]=.01
    @test op.grid!=grid
    identity=collinear_lo_operator(op.grid,2.4,2.4,profile)
    @test apply_operator(identity,h)==h
    @test apply_operator(op,h)[end]==0.
    @test dot(operator_adjoint(op,ones(length(h))),h) ≈ sum(apply_operator(op,h)) rtol=2e-15
    parameters=[.2,.3]
    f(t)=sum(apply_operator(op,(t[1]^2+exp(t[2])).*h))
    coefficient=sum(apply_operator(op,h))
    @test ForwardDiff.gradient(f,parameters) ≈ coefficient.*[2parameters[1],exp(parameters[2])] rtol=3e-15
    @test ForwardDiff.hessian(f,parameters) ≈ coefficient.*[2. 0.;0. exp(parameters[2])] rtol=3e-15
    @test_throws ErrorException ForwardDiff.derivative(q->sum(apply_operator(collinear_lo_operator(op.grid,q,2.4,profile),h)),10.)
    saved=op.matrix[1,2];op.matrix[1,2]=nextfloat(saved)
    @test_throws ErrorException apply_operator(op,h)
    @test_throws ErrorException apply_operator(op,reshape(h,:,1))
    @test_throws ErrorException operator_adjoint(op,h)
    op.matrix[1,2]=saved
    @test apply_operator(op,h)==apply_operator(old,h)
    op.grid[2]=nextfloat(op.grid[2])
    @test_throws ErrorException apply_operator(op,h)
end
