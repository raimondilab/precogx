from flask import Flask, render_template, request, jsonify
import json

app = Flask(__name__)

## Route to home page
@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    return render_template('home.html')

@app.route('/help', methods=['GET', 'POST'])
def help():
    if request.method == 'POST':
        data = request.get_json(force=True)
        print (data['gpcr'])
        return jsonify({'status': 'OK'})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Route to output page
@app.route('/output', methods=['GET', 'POST'])
def output():
    if request.method == 'POST':
        #Path to sample output
        path_to_json_output = "static/OL820/out.json"
        path_to_fasta = "static/OL820/temp/new_fasta_file.txt"
        return render_template('embedded.html', name='name', path_to_json_output=json.dumps(path_to_json_output), path_to_fasta=json.dumps(path_to_fasta))
    return ("<html><h3>It was a GET request</h3></html>")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
