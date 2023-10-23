var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TransmissionMechanisms","category":"page"},{"location":"#TransmissionMechanisms","page":"Home","title":"TransmissionMechanisms","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TransmissionMechanisms.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TransmissionMechanisms]","category":"page"},{"location":"#TransmissionMechanisms.collect_and_terms-Tuple{Any}","page":"Home","title":"TransmissionMechanisms.collect_and_terms","text":"collect_and_terms(expr::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}}; rev=true)\n\nCollect unique AND terms from a boolean expression or a vector of boolean expressions.\n\nThis function extracts unique AND terms from the input boolean expression(s) and  returns them as a vector of SymbolicUtils.BasicSymbolic{Bool}. The input expr  can be either a single boolean expression or a vector of boolean expressions.\n\nArguments\n\nexpr:  A single boolean expression or a vector of boolean expressions from which AND  terms will be collected.\nrev::Bool=true: (Optional) If true, the resulting AND terms are sorted in  descending order based on their string representation. If false, the sorting  is in ascending order.\n\nReturns\n\nReturns a vector of SymbolicUtils.BasicSymbolic{Bool} containing unique AND terms  extracted from the input expression(s).\n\nExamples\n\njulia s = \"x2 & x3 & !x2 & x5\" cond = make_condition(s) and_terms = TransmissionMechanisms.collect_and_terms(cond)`\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.contains_nots-Union{Tuple{SymbolicUtils.BasicSymbolic{T}}, Tuple{T}} where T","page":"Home","title":"TransmissionMechanisms.contains_nots","text":"contains_nots(term::SymbolicUtils.BasicSymbolic{Bool})\n\nCheck whether a term still contains a NOT statement. \n\nTo be able to calculate transmission effects using IRFs, Boolean conditions need to be simplified to a state in which they only involve ANDs in all terms. Thus, this function checks whether something went wrong in the simplification process. \n\nArguments\n\nterm::SymbolicUtils.BasicSymbolic{Bool}: A term as returned by helper_Q or get_terms. \n\nReturns\n\nReturns true if the term includes any NOT and false otherwise. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.create_transmission_function-Tuple{Int64, Int64, SymbolicUtils.BasicSymbolic{Bool}}","page":"Home","title":"TransmissionMechanisms.create_transmission_function","text":"create_transmission_function(from::Int, [to::Int,] condition::SymbolicUtils.BasicSymbolic{Bool})\n\nCreate a transmission function based on a boolean condition.\n\nThis function generates a function that calculates the transmission effect from the fromth variables to the optional toth variable. If no to is given, the effect is calculated for all destination nodes. The transmission paths must satisfy the Boolean statement given in condition which can be created using make_condition.\n\nArguments\n\nfrom::Int: Index of the variable from which the transmission starts.\n[to::Int]: (Optional) Index of the variable to which the transmission ends.  If omitted, the returned function calcuates the effect for all destination nodes. \ncondition::SymbolicUtils.BasicSymbolic{Bool}: Boolean condition specifying  the transmission mechanism. The condition is represented as a boolean expression  using variables starting with x and can be created using make_condition.\n\nReturns\n\nReturns a transmission function that can be applied to impulse response functions  (irfs) and orthogonalized impulse response functions (irfs_ortho) to compute  the transmission effects. The first argument in the returned function is irfs and the second arguement is irfs_ortho. Both should be matrices corresponding the the structural and orthogonalised IRFs respectively. Both should be square and of dimension k*(h+1) with h being the maximum horizon and k being the number of variables. The structural version can, for example, be obtained via irfs = inv(I - B) * Qbb, with B and Qbb obtained from to_structural_transmission_model.\n\nExamples\n\njulia s = \"x2 & !x3\" cond = make_condition(s) tf = create_transmission_function(1, cond) irfs = randn(3, 3) irfs_ortho = randn(3, 3) tf(irfs, irfs_ortho)`\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.get_terms-Tuple{SymbolicUtils.BasicSymbolic{Bool}}","page":"Home","title":"TransmissionMechanisms.get_terms","text":"get_terms(condition::SymbolicUtils.BasicSymbolic{Bool})\n\nSimplify a Boolean condition into a sum of only conjunctions. \n\nGiven a valid Boolean condition in the context of transmission mechanisms, the condition is recursively simplified until the result only consists of additive terms with each remaining condition being purely a conjunction. For example, The condition x2 & !x3 corresponds to the transmission statement Q(x2 & !x3) which is simplified to Q(x2) - Q(x2 & x3). The function will then return [Q(x2), -Q(x2 & x3)]. Each term can then be easily evaluated using IRFs.  \n\nArguments\n\ncondition::SymbolicUtils.BasicSymbolic{Bool}: A symbolic boolean condition. Can be created using make_condition. Valid variable names are x followed by a number, or T respresenting true. All transmission paths must satisfy this boolean statement. \n\nReturns\n\nReturns a vector of terms, each of which consists of a conjunction. \n\nExamples\n\ns = \"x2 & !x3\"\ncond = make_condition(s)\nterm = get_terms(cond)  # [Q(x2), -Q(x2 & x3)]\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.helper_Q-Tuple{SymbolicUtils.BasicSymbolic{Bool}}","page":"Home","title":"TransmissionMechanisms.helper_Q","text":"helper_Q(condition::SymbolicUtils.BasicSymbolic{Bool})\n\nSimplifies the valid Boolean condition in the context of transmission mechanisms and returns all necessary information to evaluate the transmission effect using information in Impulse Response Functions (IRFs). \n\nArguments\n\ncondition::SymbolicUtils.BasicSymbolic{Bool}: A valid symbolic Boolean condition. Can be created using make_condition. \n\nReturns\n\nA vector containing all the additive terms in the simplified representation. Each term only involves conjunctions and can thus be calculated using information in IRFs. \nA vector of multipliers applied to each term. \nA vector of variables involved in the transmission query.\nA vector of vectors. Each inner vector corresponds to a term and provides the variable numbers involved in the term. These, together with the multipliers, can be used to calculate the transmission effect. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.helper_sym_to_num-Tuple{Any}","page":"Home","title":"TransmissionMechanisms.helper_sym_to_num","text":"helper_sym_to_num(sym_vec::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}})\n\nConvert symbolic variables to numerical indices.\n\nThis function takes a single boolean symbolic variable or an array of boolean symbolic variables sym_vec and converts them into numerical indices. Each index corresponds to the variable index in the structural transmission model and thus to an index in the IRF matrix. Valid variable names should start with 'x' followed by a number or be the symbol 'T'. If invalid variable names are encountered, an error is raised.\n\nArguments\n\nsym_vec::Union{Any, Vector{Any}}:  A single boolean symbolic variable or an array of boolean symbolic variables  representing variables to be converted to numerical indices.\n\nReturns\n\nReturns a sorted vector of numerical indices corresponding to the input symbols.  If a symbol is 'T', it is represented as nothing in the output vector.\n\nExamples\n\n```julia @syms x1::Bool x2::Bool x3::Bool T::Bool y::Bool TransmissionMechanisms.helpersymtonum(x2)          # Returns [2] TransmissionMechanisms.helpersymtonum(T)          # Returns [nothing] TransmissionMechanisms.helpersymtonum([x1, x3]) # Returns [1, 3] TransmissionMechanisms.helpersymtonum(y)            # Error: Variable names should start with 'x' followed by a number.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.is_not_valid_variable_name-Tuple{SymbolicUtils.BasicSymbolic{Bool}}","page":"Home","title":"TransmissionMechanisms.is_not_valid_variable_name","text":"is_not_valid_variable_name(sym::SymbolicUtils.BasicSymbolic{Bool})\n\nCheck if a given symbol is a valid variable name in the context of transmission mechanisms.\n\nThis function validates whether the provided symbol sym is a valid variable name for use in transmission mechanisms and in create_transmission_function. In the context of transmission mechanisms, valid variable names are of the form 'x' followed by one or more digits, or the symbol 'T' representing true.\n\nArguments\n\nsym::SymbolicUtils.BasicSymbolic{Bool}: Symbol to be checked for validity as a variable name in transmission mechanisms.\n\nReturns\n\nReturns true if the input symbol is not a valid variable name, and false otherwise.\n\nExamples\n\n```julia @syms x1::Bool T::Bool y::Bool TransmissionMechanisms.isnotvalidvariablename(x1)   # Returns false (valid variable name) TransmissionMechanisms.isnotvalidvariablename(y)    # Returns true (invalid variable name) TransmissionMechanisms.isnotvalidvariablename(T)    # Returns false (valid variable name)\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_condition-Tuple{String}","page":"Home","title":"TransmissionMechanisms.make_condition","text":"make_condition(string::String)\n\nGiven a string, make a SymbolicUtils.jl Boolean expression. \n\nTransmission mechanisms are described using Boolean statements involving the  variables in the system. Each variable starts with x followed by a number. For example, given a three variable VAR(1), y{1,t} -> x1, y{2, t} -> x2,  y{3, t} -> x3, y{1, t+1} -> x4, y{2, t+1} -> x5, ... Boolean statements  then involve expressions in the x variables and define which paths can be taken. Each path involved in the transmission mechanism must satisfy the Boolean statement. \n\nArguments\n\nstring::String: A Boolean statement given as a string. Variables must start with x for them to be valid variables. \n\nReturns\n\nReturns a SymbolicUtils.jl Boolean expression that can be used in create_transmission_function. \n\nExamples\n\ns = \"x2 & !x3\"\ncond = make_condition(s)\ntf = create_transmission_function(1, cond)\nirfs = randn(3, 3)\nirfs_ortho = randn(3, 3)\ntf(irfs, irfs_ortho)\n\nNotes\n\nEnsure that variables in the input string are correctly formatted as x followed by a number. \nThe resulting Boolean expression can be used in create_transmission_function to create a function that calculates the transmission effect. \n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_structural_B-Tuple{AbstractMatrix, AbstractMatrix, Int64, Int64}","page":"Home","title":"TransmissionMechanisms.make_structural_B","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in Wegner et al (2024). \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here. \nAplus_ortho::AbstractMatrix: The structural lag matrix of the SVAR(p) obtained using a Cholesky decomposition. \np::Int: The order of the SVAR(p). \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns B in the above structural representation of the SVAR(p) over the next h horizons. \n\nNotes\n\nSee also make_structural_Qbb and to_structural_transmission_model.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.make_structural_Qbb-Tuple{AbstractMatrix, AbstractMatrix, Int64}","page":"Home","title":"TransmissionMechanisms.make_structural_Qbb","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in Wegner et al (2024). \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nirf0::AbstractMatrix: The identified impact responses. If the ith shock has not been identified, then the ith column consists of NaN. Matrix should be of dimension k times k with k being the number of variables in the SVAR(p).\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here. \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns mathbbQ in the above structural representation of the SVAR(p) over the next h horizons. \n\nNotes\n\nSee also make_structural_B and to_structural_transmission_model.\n\n\n\n\n\n","category":"method"},{"location":"#TransmissionMechanisms.to_structural_transmission_model-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Int64, Int64}","page":"Home","title":"TransmissionMechanisms.to_structural_transmission_model","text":"make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)\n\nCreate the B matrix in the structural representation in Wegner et al (2024). \n\nGiven a SVAR(p) in the form of \n\ny_tA_0 = y_t-1A_1 + dots + y_t-pA_p + varepsilon_t = y_t-1 dots y_t-pA^+ + varepsilon_t\n\nand a pre-specified maximum horizon of h, the SVAR(p) over the next h horizons can be written as \n\nx = Bx + mathbbQepsilon\n\nwith the structure of B and mathbbQ given in the paper. \n\nArguments\n\nirf0::AbstractMatrix: The identified impact responses. If the ith shock has not been identified, then the ith column consists of NaN. Matrix should be of dimension k times k with k being the number of variables in the SVAR(p).\nA0_ortho::AbstractMatrix: The contemporaneous matrix of the SVAR(p) obtained using a Cholesky decomposition. Note that this would be an upper-triangular matrix in the representation chosen here.\nAplus_ortho::AbstractMatrix: The structural lag matrix of the SVAR(p) obtained using a Cholesky decomposition. \np::Int: The order of the SVAR(p). \nmax_horizon::Int: The maximum horizon to consider for the transmission model. \n\nReturns\n\nReturns B and mathbbQ in the above structural representation of the SVAR(p) over the next h horizons in this order. \n\n\n\n\n\n","category":"method"}]
}
