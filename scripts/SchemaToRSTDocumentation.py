import os
import re
import shutil
import argparse
from xml.etree import ElementTree as etree


class TreeBuilderWithComments(etree.TreeBuilder):

    def comment(self, data):
        self.start(etree.Comment, {})
        self.data(data)
        self.end(etree.Comment)


def writeTableRST(type_name, file_name, title_prefix, values):

    element_header = '%s: %s' % (title_prefix, type_name)

    L = [[len(x) for x in row] for row in values]

    # np isn't in the docker images for our CI
    MAX = [None] * len(L[0])
    for ii in range(0, len(L[0])):
        MAX[ii] = 0
        for jj in range(0, len(L)):
            MAX[ii] = max(MAX[ii], L[jj][ii])

    M = tuple(MAX)

    # Build column formats (the final column is allowed to use multiple lines)
    single_row_format = ('{:<%is} ' * len(M) + '\n') % M

    # The final column is allowed to span multiple lines
    description_index = sum(M[:-1]) + len(M) - 1
    multiple_row_format_a = ('{:<%is} ' * (len(M) - 1) + '| {:<%is} \n') % M
    multiple_row_format_b = ' ' * description_index + '| {:<%is} \n' % (M[-1])

    # Build headers
    boundary = single_row_format.format(*['=' * x for x in M])

    # Format lines
    formatted_lines = []
    for ii in range(0, len(values)):
        if ('\\n' not in values[ii][-1]):
            # Single line entry
            formatted_lines.append(single_row_format.format(*values[ii]))

        else:
            # Multi-line entry
            description = values[ii][-1].split('\\n')
            values[ii][-1] = description[0]
            formatted_lines.append(multiple_row_format_a.format(*values[ii]))

            for d in description[1:]:
                formatted_lines.append(multiple_row_format_b.format(d))

    # Build table
    with open(file_name, 'w') as f:
        f.write('%s\n' % (element_header))
        f.write('=' * len(element_header) + '\n')
        f.write('\n')
        f.write(boundary)
        f.write(formatted_lines[0])
        f.write(boundary)
        for line in formatted_lines[1:]:
            f.write(line)
        f.write(boundary)
        f.write('\n\n')


def format_type_data(x):
    # Make sure that there are spaces within angled brackets
    # This allows for easier table line-breaks
    x = re.sub(r"<([a-zA-Z])", r"< \1", x)
    x = re.sub(r"([a-zA-Z])>", r"\1 >", x)

    # Remove any source tags
    x = re.sub(r"[a-zA-Z]*::", "", x)

    return (x)


def buildAttributeMap(root_node, xsd='{http://www.w3.org/2001/XMLSchema}'):
    attribute_map = {}

    # Build forward maps on the first pass (name, type, default, use, description, registered_by)
    for child_node in root_node.findall(xsd + 'complexType'):
        type_name = child_node.get('name')[:-4]
        attribute_map[type_name] = {}

        # Read the values supported by the xsd format (name, type, default, use)
        for attribute_node in child_node.findall(xsd + 'attribute'):
            att_name = attribute_node.get('name')
            att_type = format_type_data(attribute_node.get('type', default=' '))
            att_default = attribute_node.get('default', default=' ')

            # Handle the special case for required values
            useValue = attribute_node.get('use')
            if useValue:
                att_default = useValue

            attribute_map[type_name][att_name] = {'Type': att_type, 'Default': att_default}

        # Parse nodes
        for choice_node in child_node.findall(xsd + 'choice'):
            for sub_node in choice_node.findall(xsd + 'element'):
                sub_name = sub_node.get('name')
                sub_required = sub_node.get('minOccurs')
                sub_unique = sub_node.get('maxOccurs')
                node_use = ''
                if sub_required:
                    if sub_unique:
                        node_use = 'unique, required'
                    else:
                        node_use = 'required'
                elif (sub_unique):
                    node_use = 'unique'

                attribute_map[type_name][sub_name] = {'Type': 'node', 'Default': node_use}

        # Read the remaining values stored in comments (description, registered_by)
        for comment_node in child_node.iter(etree.Comment):
            tmp = comment_node.text.split(' => ')
            att_name = tmp[0]
            att_description = tmp[1].replace('\\\\', '\\').replace('\n', '\\n').replace('|', '\\|')
            attribute_map[type_name][att_name]['Description'] = att_description
            if (len(tmp) > 2):
                attribute_map[type_name][att_name]['Registered By'] = [x for x in tmp[2].split(', ')]

    # Build the reverse maps on the second pass (registered_on)
    for ka in attribute_map.keys():
        for kb in attribute_map[ka].keys():
            if ('Registered By' in attribute_map[ka][kb]):
                for kc in attribute_map[ka][kb]['Registered By']:
                    attribute_map[kc][kb] = attribute_map[ka][kb].copy()
                    attribute_map[kc][kb]['Registered On'] = [ka]
                    del attribute_map[kc][kb]['Registered By']

    return attribute_map


def buildTableValues(type_map, link_string='XML', include_defaults=True):
    table_headers = ['Name', 'Type', 'Default', 'Registered On', 'Registered By', 'Description']
    optional_headers = ['Default', 'Registered On', 'Registered By']

    # Sort the keys
    akeys = []
    rb_keys = []
    ro_keys = []
    node_keys = []
    for k in sorted(type_map.keys()):
        if ('Registered By' in type_map[k]):
            rb_keys.append(k)
        elif ('Registered On' in type_map[k]):
            ro_keys.append(k)
        elif (type_map[k]['Type'] == 'node'):
            node_keys.append(k)
        else:
            akeys.append(k)
    akeys.extend(rb_keys)
    akeys.extend(ro_keys)
    akeys.extend(node_keys)

    # Remove any unused headers
    for ii in range(0, len(optional_headers)):
        header_count = 0
        for k in akeys:
            if (optional_headers[ii] in type_map[k]):
                header_count += 1
        if (header_count == 0):
            table_headers.remove(optional_headers[ii])

    if not include_defaults:
        if 'Default' in table_headers:
            table_headers.remove('Default')

    # Setup the empty table
    N_rows = len(akeys)
    N_cols = len(table_headers)
    table_values = [[' ' for jj in range(0, N_cols)] for ii in range(0, N_rows + 1)]
    table_values[0] = table_headers

    if (N_rows > 0):
        # Add values to the table
        for ii in range(0, N_rows):
            # Get the row, set the name
            att_name = akeys[ii]
            table_row = table_values[ii + 1]
            table_row[0] = att_name

            # Set the other parameters
            for jj in range(1, len(table_headers)):
                k = table_headers[jj]
                if k in type_map[att_name]:
                    table_row[jj] = type_map[att_name][k]

                    # Fix type strings
                    if ('Type' in k):
                        table_row[jj] = table_row[jj].replace('_lt_', '<').replace('_gt_', '>').replace('_cm_', ',').replace('-', ' ')

                    # Format any registration entries as links
                    if ('Registered' in k):
                        table_row[jj] = ", ".join([':ref:`%s_%s`' % (link_string, x) for x in table_row[jj]])

            # Set the link if the target is a node
            if (table_row[1] == 'node'):
                table_row[-1] = ':ref:`%s_%s`' % (link_string, table_row[0])
    else:
        # Case for an empty table
        table_row = [' ' for x in table_headers]
        table_row[-1] = '(no documentation available)'
        table_values.append(table_row)

    return table_values


def main(schema_name='schema.xsd', output_folder='./', xsd='{http://www.w3.org/2001/XMLSchema}'):
    """
    Build RST Documentation Tables
    """
    # Setup folders
    additional_documentation_name = f'{schema_name}.other'
    output_folder = os.path.abspath(os.path.expanduser(output_folder))
    datastructure_folder = os.path.join(output_folder, 'datastructure')
    summary_file = os.path.join(datastructure_folder, 'CompleteXMLSchema.rst')
    os.makedirs(datastructure_folder, exist_ok=True)
    shutil.copy(schema_name, os.path.join(datastructure_folder, 'schema.xsd'))
    shutil.copy(additional_documentation_name, os.path.join(datastructure_folder, 'schema.xsd.other'))

    # Parse the input/non-input schemas
    parser = etree.XMLParser(target=TreeBuilderWithComments())
    include_tree = etree.parse(schema_name, parser=parser)
    include_root = include_tree.getroot()
    input_attribute_map = buildAttributeMap(include_root)

    parser = etree.XMLParser(target=TreeBuilderWithComments())
    include_tree = etree.parse(additional_documentation_name, parser=parser)
    include_root = include_tree.getroot()
    other_attribute_map = buildAttributeMap(include_root)

    # Check for non-unique (ignoring case) links
    input_keys = sorted(input_attribute_map.keys())
    input_keys_lower = [k.lower() for k in input_keys]
    input_keys_count = [input_keys_lower.count(k) for k in input_keys_lower]
    input_repeated_keys = [input_keys[ii] for ii in range(0, len(input_keys)) if input_keys_count[ii] > 1]

    other_keys = sorted(other_attribute_map.keys())
    other_keys_lower = [k.lower() for k in other_keys]
    other_keys_count = [other_keys_lower.count(k) for k in other_keys_lower]
    other_repeated_keys = [other_keys[ii] for ii in range(0, len(other_keys)) if other_keys_count[ii] > 1]

    if ((len(input_repeated_keys) > 0) | (len(other_repeated_keys) > 0)):
        print('Duplicate input documentation table names:')
        print(input_repeated_keys)
        print('Duplicate other documentation table names:')
        print(other_repeated_keys)
        raise ValueError('Duplicate data structure names are not allowed due to .rst limitations (case-insensitive)!')

    # Build documentation tables
    with open(summary_file, 'w') as output_handle:
        # Write the file header
        output_handle.write('###################\n')
        output_handle.write('Datastructure Index\n')
        output_handle.write('###################\n\n')

        # Parse the input schema definitions
        output_handle.write('************************\n')
        output_handle.write('Input Schema Definitions\n')
        output_handle.write('************************\n\n')

        output_handle.write(':download:`XML Schema <schema.xsd>`\n\n')

        for type_name in sorted(input_attribute_map.keys()):
            # Write the individual tables
            table_values = buildTableValues(input_attribute_map[type_name])
            writeTableRST( type_name, '%s/%s.rst' % (datastructure_folder, type_name), 'XML Element', table_values)

            # Write to the master list
            output_handle.write('\n.. _XML_%s:\n\n' % (type_name))
            output_handle.write('.. include:: %s.rst\n\n' % (type_name))

        # Parse the non-input schema definitions
        output_handle.write('*************************\n')
        output_handle.write('Datastructure Definitions\n')
        output_handle.write('*************************\n\n')

        for type_name in sorted(other_attribute_map.keys()):
            # Write the individual tables
            table_values = buildTableValues(other_attribute_map[type_name],
                                            link_string='DATASTRUCTURE',
                                            include_defaults=False)
            writeTableRST( type_name, '%s/%s_other.rst' % (datastructure_folder, type_name), 'Datastructure', table_values)

            # Write to the master list
            output_handle.write('\n.. _DATASTRUCTURE_%s:\n\n' % (type_name))
            output_handle.write('.. include:: %s_other.rst\n\n' % (type_name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--schema', type=str, help='GEOS schema file', default='schema.xsd')
    parser.add_argument('-o', '--output', type=str, help='Output directory', default='./')
    args = parser.parse_args()
    main(schema_name=args.schema, output_folder=args.output)
