#ifndef OBJECT_H
#define OBJECT_H
enum label_t
{
    SIGNAL,
    BACKGROUND
};

struct Object
{
    double* attrs;
    double weight;
    unsigned int dim;
    label_t label;
};

class ObjectCompare
{
    public:

        ObjectCompare(unsigned int attr = 0):
            attr(attr)
        {}

        bool operator()(const Object& left, const Object& right)
        {
            return left.attrs[attr] < right.attrs[attr];
        }

    private:

        unsigned int attr;
};

#endif /* OBJECT_H_ */
