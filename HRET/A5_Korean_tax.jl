# Korean tax system 

module Korean_tax

export income_tax, property_tax, acquisition_tax, HRET, rent_deduct, mortgage_lumpsum, inheritance_tax

function income_tax(earning)
    # unit = 10,000 Won
    
    # tax deduction
    if earning <= 500
        deduction = earning * 0.70
    elseif earning <= 1500
        deduction = 350 + (earning-500) * 0.40
    elseif earning <= 4500
        deduction = 750 + (earning-1500) * 0.15
    elseif earning <= 10000
        deduction = 1200 + (earning-4500) * 0.05
    else
        deduction = 1475 + (earning-10000) * 0.02
    end
    
    tax_base = earning - deduction
        
    if tax_base <= 1400
        tax = tax_base * 0.06
    elseif tax_base <= 5000
        tax = 84 + (tax_base - 1400) * 0.15
    elseif tax_base <= 8800
        tax = 624 + (tax_base - 5000) * 0.24
    elseif tax_base <= 15000
        tax = 1536 + (tax_base - 8800) * 0.35
    elseif tax_base <= 30000
        tax = 3706 + (tax_base - 15000) * 0.38
    elseif tax_base <= 50000
        tax = 9406 + (tax_base - 30000) * 0.40
    elseif tax_base <= 100000
        tax = 17406 + (tax_base - 50000) * 0.42
    else
        tax = 38406 + (tax_base - 100000) * 0.45
    end
    
    effective_tax_rate = tax / earning
    return effective_tax_rate
end

function property_tax(house_value)
    # unit = 10,000 Won
    
    # fair value adjustment
    tax_base = house_value * 0.60
    
    if tax_base <= 0
        return 0
    elseif tax_base <= 6000
        tax = tax_base * 0.001
    elseif tax_base <= 15000
        tax = 6 + (tax_base - 6000) * 0.0015
    elseif tax_base <= 30000
        tax = 19.5 + (tax_base - 15000) * 0.0025
    else
        tax = 57 + (tax_base - 30000) * 0.004
    end
    
    effective_tax_rate = tax / house_value
    return effective_tax_rate
end

function acquisition_tax(tax_base)
    # unit = 10,000 Won
    if tax_base <=0
        return 0
    elseif tax_base <= 60000
        tax = tax_base * 0.01
    elseif tax_base <= 90000
        tax = 600 + (tax_base - 60000) * 0.02
    else
        tax = 1200 + (tax_base - 90000) * 0.03
    end
    
    effective_tax_rate = tax / tax_base
    return effective_tax_rate
end

function HRET(house_value, T)
    # unit = 10,000 Won
    # tau_p = property tax
    # T     = age
    
    # deduction
    tax_base = (house_value-60000)*0.60
    
    if tax_base <=0
        tax = 0
        return 0
    elseif tax_base <=60000
        tax = tax_base*0.005
    elseif tax_base <= 120000
        tax = 300 + (tax_base-60000)*0.0075
    elseif tax_base <= 500000
        tax = 750 + (tax_base-120000)*0.01
    elseif tax_base <= 940000
        tax = 4550 + (tax_base-500000)*0.015
    else
        tax = 11150 + (tax_base-940000)*0.02
    end
    
    # property tax
    tau_p_common = (house_value-60000)*0.60*0.004
    tax_p_adj = tax - tau_p_common
    
    # final tax
    if tax_p_adj<=0
        final_tax = 0
        return 0
    elseif T < 36
        final_tax = tax_p_adj
    elseif T < 41 # (60~64)
        final_tax = tax_p_adj*0.80
    elseif T < 46 # (65~69)
        final_tax = tax_p_adj*0.70
    else
        final_tax = tax_p_adj*0.60
    end
    
    effective_tax_rate = final_tax / house_value
    return effective_tax_rate
    
end

function rent_deduct(earning)
    # unit = 10,000 Won
    # rent = ρ*s
    if earning <=5500
        return 0.17
    elseif earning <= 7000
        return 0.15
    else
        return 0
    end
end

function mortgage_lumpsum(house_value, debt, r_m)
    # unit = 10,000 Won
    # debt is negative value
    # rent = ρ*s
    if house_value <=50000
        return -r_m*debt
    else
        return 0
    end
end

function inheritance_tax(tax_base)
    if tax_base <=0
        tax = 0
        return 0
    elseif tax_base <=10000
        tax = tax_base*0.10
    elseif tax_base <=50000
        tax = 1000 + (tax_base-10000)*0.20
    elseif tax_base <= 100000
        tax = 9000 + (tax_base-50000)*0.30
    elseif tax_base <= 300000
        tax = 24000 + (tax_base-100000)*0.40
    else
        tax = 104000 + (tax_base-300000)*0.50
    end
    avg_tax_rate = tax / tax_base
    return avg_tax_rate
end

end