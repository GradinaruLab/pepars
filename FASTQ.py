
def convert_quality_scores_to_probability_of_error(quality_scores):
	
	probability_of_error_vector = []

	for quality_score in quality_scores:
		probability_of_error = convert_quality_score_to_probability_of_error(quality_score)
		probability_of_error_vector.append(probability_of_error)

	return probability_of_error_vector

def convert_quality_score_to_probability_of_error(quality_score):
	
	return pow(10,quality_score/(-10))
	
def convert_quality_string_to_quality_score(quality_score_string):
	quality_score_vector = []
	
	for current_character in quality_score_string:
		#Convert character into a number and subtract 33 to compensate for illumina ascii format
		quality_score = ord(current_character)-33
		quality_score_vector.append(quality_score)
		
	return quality_score_vector