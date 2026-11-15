"An explicitly supplied qT=P_hperp/z1 window at fixed theory z1."
struct SIAQTWindow
    lower_GeV::Float64
    upper_GeV::Float64
    source::String
    lower_limit_source::String
    function SIAQTWindow(lower::Float64,upper::Float64,source::String,lower_source::String)
        isfinite(lower) && isfinite(upper) && 0<=lower<upper ||
            error("finite explicit qT bounds 0<=lower<upper required")
        !isempty(strip(source)) && !isempty(strip(lower_source)) ||
            error("nonblank qT-window source and lower-limit provenance required")
        new(lower,upper,source,lower_source)
    end
end
SIAQTWindow(lower,upper;source,lower_limit_source)=SIAQTWindow(
    fixed_scalar(lower,"fixed lower qT"),fixed_scalar(upper,"fixed upper qT"),
    String(source),String(lower_limit_source))

"""
Prepare the parent's separate P*dP yields over P=z1*qT. Both endpoints and
their provenance are caller supplied. No event/bin acceptance or experimental
admission follows from the coordinate conversion or a reported mean z1.
"""
function sia_qt_observable(context::PolarizedLeptonicPlan,window::SIAQTWindow;
        ratio,angular::Union{SIAAngularFactor,SIAMeanOfRatio},order=32)
    context.leptonic.process==:SIA || error("SIA hadron2 context required for qT window")
    z1=context.leptonic.z1
    lower=z1*window.lower_GeV;upper=z1*window.upper_GeV
    window.lower_GeV>0 && lower==0 && error("mapped positive lower P bound underflows")
    momentum=MomentumRule(lower,upper;order=order)
    metadata=Dict{String,Any}(
        "frame"=>"hadron2","double_ratio_approximation"=>"linearized",
        "qT_window_contract_id"=>"KPSY-SIA-QT-WINDOW-001",
        "input_momentum_coordinate"=>"qT=P_hperp/z1 at fixed theory z1",
        "input_momentum_bounds_GeV"=>[window.lower_GeV,window.upper_GeV],
        "input_momentum_measure"=>"z1^2*qT*dqT = P_hperp*dP_hperp",
        "fixed_theory_z1"=>z1,"qT_window_source"=>window.source,
        "qT_lower_limit_source"=>window.lower_limit_source,
        "kinematics"=>"fixed_other_kinematics_finite_qT_integral",
        "acceptance"=>"uniform momentum acceptance; explicit reconstruction assumption only",
        "experimental_data_admitted"=>false,"historical_fit_ready"=>false,"G2PERF_pass"=>false,
        "joint_acceptance_recovered"=>false,"energy_momentum_fraction_conversion"=>"not supplied")
    merge!(metadata,sia_angular_metadata(angular))
    assemble_observable(ratio,context,momentum,angular.value,window.source;extra=metadata)
end
