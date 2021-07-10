# Qualification - Problem 2.2
This is a relatively simple problem once understood (in my opinion, it was confusingly worded).
With each sequencing treated as a "row", the problem is basically counting the total number of unique "columns" and assign them an ID.
With that concept in mind, I basically joined the column into a string key that mapped to the next numerical ID in a dictionary.
Any new keys were given assigned the next available ID.
