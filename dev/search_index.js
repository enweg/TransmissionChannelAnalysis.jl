var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TransmissionMechanisms","category":"page"},{"location":"#TransmissionMechanisms","page":"Home","title":"TransmissionMechanisms","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TransmissionMechanisms.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TransmissionMechanisms]","category":"page"},{"location":"#TransmissionMechanisms.Q","page":"Home","title":"TransmissionMechanisms.Q","text":"Q(s::String)\nQ(s::String, m::Number)\nQ(s::Vector{String})\nQ(s::Vector{String}, m::Vector{Number})\nQ(i::Int)\n\nRepresents a transmission condition. \n\nWe denote with Q(b), where b is a Boolean statement, a transmission question. \n\nImportant: Each string in s: should be a Boolean statement involving variables xfollowed by a number, i.e.x1,x2`, etc. Additionally, each Boolean statement should contain only AND (&) and NOT (!) statements. \n\nImportant: Users should only use the Q(i::Int) constructor. All other constructors are for internal use only. Misuse of the other constructors easily leads to mistakes. \n\nFields\n\nvars::Vector{String}: Contains the variables. These will be Boolean statements containing only AND and NOT.\nmultiplier::Vector{Number}: Multiplier in front of Q(b). \n\nArguments\n\ns::Union{String, Vector{String}}: String representation of transmission condition. \nm::Union{Number, Vector{Number}}: Multipliers for transmission conditions. \ni::Int: Variable number.\n\nUsage\n\n\n# Defining all variables in one go\nx = [Q(\"x$i\") for i = 1:10]\nq = (x[1] | x[2]) & !x[3]\n\n# Alternatively variables can be defined separaterly\nx1 = Q(\"x1\")\nx2 = Q(\"x2\")\nx3 = Q(\"x3\")\nq = (x1 | x2) & !x3\n\n# The following are also valid\nq = Q(\"x1 & !x3\", 1)\nq = Q([\"x1\", \"x2\", \"x1 & x2\"], [1, 1, -1])\n\n# The following is NOT valid but does not yet throw an error or warning\nq = Q(\"x1 | x2\")  # DO NOT DO THIS!\n\n\n\n\n\n","category":"type"},{"location":"#Base.:!-Tuple{TransmissionMechanisms.Q}","page":"Home","title":"Base.:!","text":"Base.!(q1::Q)\n\nReturn NOT the transmission condition if the condition involves more than one variables. If the condition only involves one variables, then \"!x1\" is returned where \"1\" is replaced by the respective variable number. \n\nNote: The decision not to simplify terms of the form \"!x1\" was made because calculations usign the second calculation method in Wegner et al (2024) are faster than having to simplify \"!x1\" type of terms and using the first calculation method. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:&-Tuple{TransmissionMechanisms.Q, TransmissionMechanisms.Q}","page":"Home","title":"Base.:&","text":"Base.&(q1::Q, q2::Q)\n\nCombine two transmission conditions using AND. \n\n\n\n\n\n","category":"method"},{"location":"#Base.:|-Tuple{TransmissionMechanisms.Q, TransmissionMechanisms.Q}","page":"Home","title":"Base.:|","text":"Base.|(q1::Q, q2::Q)\n\nCombine two transmission conditions using OR. \n\n\n\n\n\n","category":"method"},{"location":"#Base.string-Tuple{TransmissionMechanisms.Q}","page":"Home","title":"Base.string","text":"Base.string(q::Q)\n\nObtain a string representation of a transmission condition.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.apply_and!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, AbstractMatrix{T}, Int64, Int64}} where T","page":"Home","title":"TransmissionMechanisms.apply_and!","text":"apply_and!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int)\n\nManipulate B and Qbb so that var lies on all paths. This corresponds to zeroing out all edges going directly from the shock to any variables ordered after var and zeroing out any edges going from variables ordered before var to any variables ordered after var.  \n\nArguments\n\nB::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_B. \nQbb::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_Qbb. \nfrom::Int: The shock number. \nvar::Int: The variable number that must lie on all paths. \n\nNotes\n\nThis function is meant for internal use only. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.apply_not!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, AbstractMatrix{T}, Int64, Int64}} where T","page":"Home","title":"TransmissionMechanisms.apply_not!","text":"apply_not!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int)\n\nManipulate B and Qbb so that var lies on no paths. This corresponds to zeroing out the edge from the shock to var, and zeroing out all edges from variables ordered before var to var. The paper mentions also zeroing out edges leaving var, but this is not necessary.  \n\nArguments\n\nB::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_B. \nQbb::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_Qbb. \nfrom::Int: The shock number. \nvar::Int: The variable number that must lie on all paths. \n\nNotes\n\nThis function is meant for internal use only. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.calculate_Q_and_only-Union{Tuple{T}, Tuple{Int64, AbstractMatrix{T}, AbstractMatrix{T}, AbstractVector{Int64}, Number}} where T","page":"Home","title":"TransmissionMechanisms.calculate_Q_and_only","text":"calculate_Q_and_only(from::Int, \n    irfs::AbstractMatrix{T}, \n    irfs_ortho::AbstractMatrix{T}, \n    vars::AbstractVector{Int}, \n    multiplier::Number\n) where {T}\n\nCalculate the transmission effect of a transmission condition/query that involves only ANDs. \n\nArguments\n\nfrom::Int: ID of shock. \nirfs::AbstractMatrix{T}: IRFs in transmission form. See to_transmission_irfs.\nirfs_ortho::AbstractMatrix{T}: Orthogonalised IRFs in transmission form. See to_transmission_irfs.\nmultiplier::Number: Multiplier.\n\nReturns\n\nReturns a Vector{T} with entry i corresponding to the transmission effect on variable x_i. If x_k is the variable in the transmission condition with the highest subscript, then all entries in the returned vector with index less thatn k are NaN since interpretation of those results is nonsensical.\n\nNotes\n\nOnly for internal use. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.check_contradiction-Tuple{Vector{Int64}, Vector{Int64}}","page":"Home","title":"TransmissionMechanisms.check_contradiction","text":"check_contradiction(var_and::Vector{Int}, var_not::Vector{Int})\ncheck_contradiction(var_and::Vector{Vector{Int}}, var_not::Vector{Vector{Int}})\n\nCheck whether there is a contradiction of the form x1 & !x1. \n\nArguments\n\nvar_and::Union{Vector{Int}, Vector{Vector{Int}}}: AND variable numbers obtained from get_varnums_and_multiplier. \nvar_not::Union{Vector{Int}, Vector{Vector{Int}}}: NOT variable numbers obtained from get_varnums_and_multiplier\n\nReturns\n\nBool indicating whether there are any contradictions. \nVector{Bool} indicating which elements yielded a contradiction.\n\nNotes\n\nThis is used in remove_contradictions to remove contradicting terms. This speeds up simplification of terms, since the total number of terms can often be reduced. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.collect_terms-Tuple{TransmissionMechanisms.Q}","page":"Home","title":"TransmissionMechanisms.collect_terms","text":"collect_terms(q::Q)\n\nCollect all terms Q(b) for which the Boolean statement b is the same and sums their multiplier. The result is a transmission condition for which each term only appears ones, but with multipliers possibly different from plus-minus one. \n\nArguments\n\nq::Q: A transmission condition. See also Q\n\nReturns\n\nAnother transmission condition of type Q. \n\nExample\n\nq = Q([\"x1\", \"x1\"], [1, 1])  \ncollect_terms(q)\n# output: Q(\"x1\", 2)\nq = Q([\"x1\", \"\", \"x1\"], [1, 1, -1])  \ncollect_terms(q)\n# output: Q(\"\", 1)\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.get_varnums_and_multiplier-Tuple{TransmissionMechanisms.Q}","page":"Home","title":"TransmissionMechanisms.get_varnums_and_multiplier","text":"get_varnums_and_multiplier(q::Q)\n\nObtain the AND and NOT expressions of a transmission condition. Also return the multiplier for each term. \n\nA valid transmission condition is a set of terms involving only AND and NOT expressions. For each term, the AND and NOT expressions are collected and a vector of the respective variable numbers is returned. \n\nArguments\n\nq::Q: A transmission condition. See also Q. \n\nReturns\n\nand_nums::Vector{Vector{Int}}: Contains for each term a vector of variable numbers that are included via AND in the term. \nand_not_nums::Vector{Vector{Int}}: Contains for each term a vector of variable numbers that are included via NOT in the term. \nmultiplier::Vector{Number}: Contains for each term the multiplier. \n\nExamples\n\nq = Q([\"x1\", \"!x2\", \"x1 & !x2\"], [1, 2, 3])\nand_nums, and_not_nums, multiplier = get_varnums_and_multiplier(q)\n# output: \n# and_nums = [[1], [], [1]]\n# and_not_nums = [[], [2], [2]]\n# multiplier = [1, 2, 3]\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_condition-Tuple{String}","page":"Home","title":"TransmissionMechanisms.make_condition","text":"make_condition(s::String)\n\nGiven a s::String, form a transmission condition of the form Q(b) where b is a Boolean statement. \n\nTransmission mechanisms are described using Boolean statements involving the  variables in the system. Each variable starts with x followed by a number. For example, given a three variable VAR(1), y{1,t} -> x1, y{2, t} -> x2,  y{3, t} -> x3, y{1, t+1} -> x4, y{2, t+1} -> x5, ... Boolean statements  then involve expressions in the x variables and define which paths can be taken. Each path involved in the transmission mechanism must satisfy the Boolean statement. \n\nArguments\n\ns::String: A Boolean statement given as a string. Variables must start with x for them to be valid variables. \n\nReturns\n\nReturns a transmission condition. See also Q.\n\nExamples\n\ns = \"x2 & !x3\"\ncond = make_condition(s)\n\nNotes\n\nEnsure that variables in the input string are correctly formatted as x followed by a number. \nThe resulting transmission condition can be used in transmission to calculate the transmission effect.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_structural_B-Tuple{AbstractMatrix, AbstractMatrix, Int64, Int64}","page":"Home","title":"TransmissionMechanisms.make_structural_B","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in WEGNER. \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here. \nAplus_ortho::AbstractMatrix: The structural lag matrix of the SVAR(p) obtained using a Cholesky decomposition. \np::Int: The order of the SVAR(p). \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns B in the above structural representation of the SVAR(p) over the next h horizons. \n\nNotes\n\nSee also make_structural_Qbb and to_structural_transmission_model.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_structural_Qbb-Tuple{AbstractMatrix, AbstractMatrix, Int64}","page":"Home","title":"TransmissionMechanisms.make_structural_Qbb","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in WEGNER. \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nirf0::AbstractMatrix: The identified impact responses. If the ith shock has not been identified, then the ith column consists of NaN. Matrix should be of dimension k times k with k being the number of variables in the SVAR(p).\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here. \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns mathbbQ in the above structural representation of the SVAR(p) over the next h horizons. \n\nNotes\n\nSee also make_structural_B and to_structural_transmission_model.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.remove_contradictions-Tuple{TransmissionMechanisms.Q}","page":"Home","title":"TransmissionMechanisms.remove_contradictions","text":"remove_contradictions(q::Q)\n\nRemove contradicting terms. \n\nA terms is deemed contradicting if it includes some \"xi & !xi\". This would result in the entire Boolean statement to be false, and thus in the effect of this terms to be zero. \n\nArguments\n\nq::Q: A transmission condition. See also Q and make_condition. \n\nReturns\n\nIf TransmissionMechanisms.REMOVE_CONTRADICTIONS == false, then q will simply be returned again. \nIf TransmissionMechanisms.REMOVE_CONTRADICTIONS == false, then \nIf all terms are contradicting, the Q(\"\", 0) will be retuned, which has a transmission effect of zero. \nIf some terms are non-contradicting, then a transmission condition consisting of only the non-contradicting terms will be returned. \n\nExamples\n\nTransmissionMechanisms.REMOVE_CONTRADICTIONS = true\nq = TransmissionMechanisms.Q(\"x1\", 1)\nremove_contradictions(q)  # will return q again since no contradictions exist\n\nq = TransmissionMechanisms.Q(\"x1 & !x1\", 1)\nremove_contradictions(q)  # Will return Q(\"\", 0)\n\nq = TransmissionMechanisms.Q([\"x1 & !x1\", \"x1 & x2\"], [1, 1])\nremove_contradictions(q)  # Will return Q(\"x1 & x2\", 1)\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.string_and-Tuple{String, String}","page":"Home","title":"TransmissionMechanisms.string_and","text":"string_and(s1::String, s2::String)\n\nCombine two strings using \"&\".\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.to_structural_transmission_model-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Int64, Int64}","page":"Home","title":"TransmissionMechanisms.to_structural_transmission_model","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in WEGNER. \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nirf0::AbstractMatrix: The identified impact responses. If the ith shock has not been identified, then the ith column consists of NaN. Matrix should be of dimension k times k with k being the number of variables in the SVAR(p).\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here.\nAplus_ortho::AbstractMatrix: The structural lag matrix of the SVAR(p) obtained using a Cholesky decomposition. \np::Int: The order of the SVAR(p). \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns B and mathbbQ in the above structural representation of the SVAR(p) over the next h horizons in this order. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.to_transmission_irfs-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T","page":"Home","title":"TransmissionMechanisms.to_transmission_irfs","text":"to_transmission_irfs(irfs::AbstractArray{T, 3})\n\nTransform a standard three dimensional IRF array into a IRF matrix. The \n\nArguments\n\nirfs::AbstractArray{T, 3}: IRF Array of dimension nvariabels × nshocks × n_horizons with the first horizons corresponding to horizon 0. \n\nReturns\n\nMatrix{T} of dimension (nvariables * nhorizons) × (nvariables * nhorizons). This is the same as what would be obtained via (I-B)^-1mathbbQ using the notation of Wegner et al (2024). \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.transmission-Union{Tuple{T}, Tuple{Int64, AbstractMatrix{T}, AbstractMatrix{T}, TransmissionMechanisms.Q, Type{Val{:BQbb}}}} where T","page":"Home","title":"TransmissionMechanisms.transmission","text":"transmission(from::Int, \n    B::AbstractMatrix{T},\n    Qbb::AbstractMatrix{T}, \n    q::Q, \n    ::Type{Val{:BQbb}}\n) where {T}\n\nGiven a transmission condition q, calculate the transmission effect using the :BQbb method. \n\nArguments\n\nfrom::Int: Shock number. \nB::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_B. \nQbb::AbstractMatrix{T}: Part of the structural transmission representation in Wegner et al (2024). See also make_structural_Qbb. \nq::Q: A transmission condition. See also Q. \n\nReturns\n\nReturns a Vector{T} with entry i corresponding to the transmission effect on variable x_i. If x_k is the variable in the transmission condition with the highest subscript, then all entries in the returned vector with index less thatn k are NaN since interpretation of those results is nonsensical. \n\nExamples\n\nk = 6\nh = 3\ns = \"(x1 | x2) & !x3\"\ncond = make_condition(s)\n\nB = randn(k*(h+1), k*(h+1))\nQbb = randn(k*(h+1), k*(h+1))\n\neffect = transmission(1, B, Qbb, cond)\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.transmission-Union{Tuple{T}, Tuple{Int64, AbstractMatrix{T}, AbstractMatrix{T}, TransmissionMechanisms.Q, Type{Val{:irfs}}}} where T","page":"Home","title":"TransmissionMechanisms.transmission","text":"transmission(from::Int, \n    irfs::AbstractMatrix{T}, \n    irfs_ortho::AbstractMatrix{T}, \n    q::Q, \n    ::Type{Val{:irfs}}\n) where {T}\n\nGiven a transmission condition q, calculate the transmission effect using the :irfs method. \n\nArguments\n\nfrom::Int: Shock number. \nirfs::AbstractMatrix{T}: Impulse responses. These should be in the form of the structural transmission model. See also to_transmission_irfs. \nirfs_ortho::AbstractMatrix{T}: Orthogonalised IRFs. These should be in the form of the structural transmission model. See also to_transmission_irfs. \nq::Q: A transmission condition. See also Q. \n\nReturns\n\nReturns a Vector{T} with entry i corresponding to the transmission effect on variable x_i. If x_k is the variable in the transmission condition with the highest subscript, then all entries in the returned vector with index less thatn k are NaN since interpretation of those results is nonsensical. \n\nExamples\n\nk = 6\nh = 3\ns = \"(x1 | x2) & !x3\"\ncond = make_condition(s)\n\nirfs = randn(k, k, h+1)\nirfs_ortho = randn(k, k, h+1)\n\nirfs = to_transmission_irfs(irfs)\nirfs_ortho = to_transmission_irfs(irfs_ortho)\n\neffect = transmission(1, irfs, irfs_ortho, cond; method = :irfs)\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.transmission-Union{Tuple{T}, Tuple{Int64, AbstractMatrix{T}, AbstractMatrix{T}, TransmissionMechanisms.Q}} where T","page":"Home","title":"TransmissionMechanisms.transmission","text":"transmission(from::Int, arr1::AbstractMatrix{T}, arr2::AbstractMatrix{T}, q::Q; method = :BQbb) where {T}\n\nGiven a transmission condition q, calculate the transmission effect using the either the :BQbb method (the default), or the :irfs method. \n\nArguments\n\nfrom::Int: Shock number. \narr1::AbstractMatrix{T}. In case of :BQbb this must be B, in case of :irfs this must be irfs. See the documentation for the specific methods for transmission(..., ::Type{Val{:BQbb}}) and transmission(...,::Type{Val{:irfs}}). \narr2::AbstractMatrix{T}: In case of :BQbb this must be Qbb, in case of :irfs this must be irfs_ortho. See the documentation for the specific methods for transmission(..., ::Type{Val{:BQbb}}) and transmission(...,::Type{Val{:irfs}}).  \nq::Q: A transmission condition. See also Q.\n\nKeyword Arguments\n\nmethod::Symbol: Either :BQbb in which case the transmission effect will be calculated using the second method in Wegner et al (2024), or :irfs in which case the transmission effect is calculated using the first method in Wegner et al (2024). \n\nReturns\n\nReturns a Vector{T} with entry i corresponding to the transmission effect on variable x_i. If x_k is the variable in the transmission condition with the highest subscript, then all entries in the returned vector with index less thatn k are NaN since interpretation of those results is nonsensical. \n\nExamples\n\nk = 6\nh = 3\ns = \"(x1 | x2) & !x3\"\ncond = make_condition(s)\n\nB = randn(k*(h+1), k*(h+1))\nQbb = randn(k*(h+1), k*(h+1))\n\neffect = transmission(1, B, Qbb, cond; method = :BQbb)\neffect = transmission(1, B, Qbb, cond)  # same as above; default is :BQbb\n\nirfs = randn(k, k, h+1)\nirfs_ortho = randn(k, k, h+1)\n\nirfs = to_transmission_irfs(irfs)\nirfs_ortho = to_transmission_irfs(irfs_ortho)\n\neffect = transmission(1, irfs, irfs_ortho, cond; method = :irfs)\n\n\n\n\n\n","category":"method"}]
}
