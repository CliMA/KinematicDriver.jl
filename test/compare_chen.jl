using CSV, DataFrames

chen_data = CSV.read("q_rai_plt_chen.csv", DataFrame, header= false)
original_data = CSV.read("q_rai_plt_original.csv", DataFrame, header = false)

print(typeof(chen_data))
print(typeof(original_data))

if Matrix(chen_data) == Matrix(original_data)
    print("equal values")
end