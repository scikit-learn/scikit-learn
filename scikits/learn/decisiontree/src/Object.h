#ifndef OBJECT_H
#define OBJECT_H

struct Object
{
    double* attrs;
    unsigned int dim;
    double weight;
    bool label;
};

class ObjectCompare
{
    public:

        ObjectCompare(unsigned int attr = 0):
            attr(attr)
        {}

        bool operator()(const Object& left, const Object& right)
        {
            return *left.attrs[attr] < *right.attrs[attr];
        }

    private:

        int attr;
};

#endif /* OBJECT_H_ */
