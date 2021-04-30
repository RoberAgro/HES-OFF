from hes_off_app import app
import webbrowser
import time

if __name__ == "__main__":
    # webbrowser.open("http://127.0.0.1:8080/", new=0)
    app.run(host="127.0.0.1", port=8080, debug=True)





# TODO
#  Get a better understanding of the routes and methods (DONE Flask)
#  Create a button that computes the performance when clicked (DONE)
#  The fields should not be reset after clicking the compute button (DONE WTForms)
#  The form should use numeric inputs with checks for min/max value (DONE WTForms)
#  I have to deal with the units of the input variables in a nicer way (DONE)
#  There is a HTML <input> to load a file. It could be nice to load the configuration file from local (DONE)
#  Buttom to grey-out and do not use WT, FC or EL and simply set zero ro rated powers instead (IGNORED)
#  Use a fieldset with a legend to group up the items related to process, GT, WT, FC, EL, H2 (IGNORED)
#  Add the required option to all the input fields (DONE WTForms)
#  Clean the problems with FC/EL_EFFICIENCY initialization (DONE)
#  Clean the main function route to make it readable (DONE)
#  Think whether to use multiple forms or have a single massive one as I am doing now (SINGLE FORM DONE)
#  Think about the structure for CSS (Nice Jinja loops + Bootrap coat of paint for formatting)
#  Add a dropdown menu with multiple selectable fields to choose what plots to show
#  Add integer validator to GT (Done)



#  I can make two columns for input section and output section
#  I can make inline blocks with different boxes to contain the label and the items
#  I have to find a nice template to style the webpage
# I have to use short variable names with a given number of characters so that it occupies half a page
# I have to learn how to display a message on hover

#  cONTINUE WITH css TUTORIALS AND LEARN HOW to make a two-column table

# Miguel tutorials on bootstrap and hovering
#  Add placeholdrs peak mid tail dynamically in html/jinja?

# TODO I have to learn more about bootstrap and create a webpage with 3 pages


# At the start of routes create two dictionaries with the keys (label,data) for figs and data
# create_data_dict()
# fill_data_dict() that uses data from evaluate model to fill in the dictionary