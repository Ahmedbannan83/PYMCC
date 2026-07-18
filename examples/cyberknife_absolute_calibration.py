from pymcc.cyberknife_calibration import calculate_from_mcc, result_as_dict

# Per Accuray Chapter 2, acquire the PDD with the 60 mm cone at SSD=1000 mm.
# The supplied PDD.mcc reports SSD=800 mm, so the result deliberately carries
# a non-compliance warning and must not be used for clinical calibration.
result = calculate_from_mcc(
    "PDD.mcc",
    ionization_chamber="PTW 31021 Semiflex 3D",
)
print(result_as_dict(result))
