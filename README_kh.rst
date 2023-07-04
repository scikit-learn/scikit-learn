scikit-learn គឺជាម៉ូឌុល Python សម្រាប់ការរៀនម៉ាស៊ីនដែលបង្កើតឡើងនៅលើកំពូលនៃ SciPy ហើយត្រូវបានចែកចាយក្រោមអាជ្ញាប័ណ្ណ 3-Clause BSD ។

គម្រោងនេះត្រូវបានចាប់ផ្តើមនៅឆ្នាំ 2007 ដោយលោក David Cournapeau ជាគម្រោង Google Summer of Code ហើយចាប់តាំងពីពេលនោះមក អ្នកស្ម័គ្រចិត្តជាច្រើនបានចូលរួមចំណែក។ សូមមើលទំព័រអំពីពួកយើង សម្រាប់បញ្ជីអ្នករួមចំណែកស្នូល។

បច្ចុប្បន្ន វាត្រូវបានថែរក្សាដោយក្រុមអ្នកស្ម័គ្រចិត្ត។

គេហទំព័រ៖ https://scikit-learn.org

ការដំឡើង
ភាពអាស្រ័យ
scikit-learn ទាមទារ៖

Python (>= 3.8)
NumPy (>= 1.17.3)
SciPy (>= 1.5.0)
joblib (>= 1.1.1)
threadpoolctl (>= 2.0.0)
Scikit-learn 0.20 គឺជាកំណែចុងក្រោយដើម្បីគាំទ្រ Python 2.7 និង Python 3.4 ។ scikit-learn 1.0 ហើយក្រោយមកទាមទារ Python 3.7 ឬថ្មីជាងនេះ។ scikit-learn 1.1 ហើយក្រោយមកទាមទារ Python 3.8 ឬថ្មីជាងនេះ។

សមត្ថភាពធ្វើផែនការ Scikit-learn (ឧ. មុខងារចាប់ផ្តើមដោយ plot_ ហើយថ្នាក់បញ្ចប់ដោយ "Display") ទាមទារ Matplotlib (>= 3.1.3)។ សម្រាប់ការដំណើរការឧទាហរណ៍ Matplotlib >= 3.1.3 ត្រូវបានទាមទារ។ ឧទាហរណ៍មួយចំនួនទាមទារ scikit-image >= 0.16.2 ឧទាហរណ៍មួយចំនួនត្រូវការខ្លាឃ្មុំផេនដា >= 1.0.5 ឧទាហរណ៍ខ្លះទាមទារ seaborn >= 0.9.0 និង plotly >= 5.14.0 ។

ការដំឡើងអ្នកប្រើប្រាស់
ប្រសិនបើអ្នកមានការដំឡើង numpy និង scipy រួចហើយ វិធីងាយស្រួលបំផុតក្នុងការដំឡើង scikit-learn គឺប្រើ pip៖

pip ដំឡើង -U scikit-learn
ឬ conda:

conda ដំឡើង -c conda-forge scikit-learn
ឯកសាររួមមានការណែនាំលម្អិតបន្ថែមអំពីការដំឡើង។

កំណត់ហេតុផ្លាស់ប្តូរ
សូមមើលកំណត់ហេតុនៃការផ្លាស់ប្តូរសម្រាប់ប្រវត្តិនៃការផ្លាស់ប្តូរគួរឱ្យកត់សម្គាល់ចំពោះ scikit-learn ។

ការអភិវឌ្ឍន៍
យើងស្វាគមន៍អ្នករួមចំណែកថ្មីនៃកម្រិតបទពិសោធន៍ទាំងអស់។ គោលដៅសហគមន៍ scikit-learn គឺមានប្រយោជន៍ ស្វាគមន៍ និងមានប្រសិទ្ធភាព។ មគ្គុទ្ទេសក៍អភិវឌ្ឍន៍មានព័ត៌មានលំអិតអំពីការរួមចំណែក លេខកូដ ឯកសារ ការធ្វើតេស្ត និងច្រើនទៀត។ យើងបានបញ្ចូលព័ត៌មានមូលដ្ឋានមួយចំនួននៅក្នុង README នេះ។

តំណភ្ជាប់សំខាន់ៗ
កូដប្រភពផ្លូវការ៖ https://github.com/scikit-learn/scikit-learn
ទាញយកការចេញផ្សាយ៖ https://pypi.org/project/scikit-learn/
កម្មវិធីតាមដានបញ្ហា៖ https://github.com/scikit-learn/scikit-learn/issues
ប្រភព​កូដ
អ្នកអាចពិនិត្យមើលប្រភពចុងក្រោយបំផុតដោយប្រើពាក្យបញ្ជា៖

ក្លូន git https://github.com/scikit-learn/scikit-learn.git
ការរួមចំណែក
ដើម្បីស្វែងយល់បន្ថែមអំពីការរួមចំណែកដល់ scikit-learn សូមមើលការណែនាំអំពីការរួមចំណែករបស់យើង។

ការធ្វើតេស្ត
បន្ទាប់ពីដំឡើងរួច អ្នកអាចបើកដំណើរការឈុតសាកល្បងពីខាងក្រៅថតប្រភព (អ្នកនឹងត្រូវដំឡើង pytest >= 7.1.2)៖

pytest sklearn
សូមមើលគេហទំព័រ https://scikit-learn.org/dev/developers/contributing.html#testing-and-improving-test-coverage សម្រាប់ព័ត៌មានបន្ថែម។

ការបង្កើតលេខចៃដន្យអាចត្រូវបានគ្រប់គ្រងកំឡុងពេលធ្វើតេស្តដោយកំណត់អថេរបរិស្ថាន SKLEARN_SEED ។
ការដាក់ស្នើរសុំទាញ
មុនពេលបើកសំណើទាញ សូមក្រឡេកមើលទំព័រការរួមចំណែកពេញលេញ ដើម្បីប្រាកដថាលេខកូដរបស់អ្នកអនុលោមតាមគោលការណ៍ណែនាំរបស់យើង៖ https://scikit-learn.org/stable/developers/index.html

ប្រវត្តិគម្រោង
គម្រោងនេះត្រូវបានចាប់ផ្តើមនៅឆ្នាំ 2007 ដោយលោក David Cournapeau ជាគម្រោង Google Summer of Code ហើយចាប់តាំងពីពេលនោះមក អ្នកស្ម័គ្រចិត្តជាច្រើនបានចូលរួមចំណែក។ សូមមើលទំព័រអំពីពួកយើង សម្រាប់បញ្ជីអ្នករួមចំណែកស្នូល។

គម្រោងនេះបច្ចុប្បន្នត្រូវបានថែរក្សាដោយក្រុមអ្នកស្ម័គ្រចិត្ត។

ចំណាំ៖ scikit-learn ពីមុនត្រូវបានគេហៅថា scikits.learn ។

ជួយ​និង​គាំទ្រ
ឯកសារ
ឯកសារ HTML (ការចេញផ្សាយប្រកបដោយស្ថេរភាព): https://scikit-learn.org
ឯកសារ HTML (កំណែអភិវឌ្ឍន៍): https://scikit-learn.org/dev/
សំណួរគេសួរញឹកញាប់៖ https://scikit-learn.org/stable/faq.html
ការ​ទំនាក់ទំនង
បញ្ជីសំបុត្ររួម៖ https://mail.python.org/mailman/listinfo/scikit-learn
Gitter៖ https://gitter.im/scikit-learn/scikit-learn
ឡូហ្គោ និងម៉ាកយីហោ៖ https://github.com/scikit-learn/scikit-learn/tree/main/doc/logos
ប្លុក៖ https://blog.scikit-learn.org
ប្រតិទិន៖ https://blog.scikit-learn.org/calendar/
Twitter៖ https://twitter.com/scikit_learn
Stack Overflow៖ https://stackoverflow.com/questions/tagged/scikit-learn
ការពិភាក្សា Github៖ https://github.com/scikit-learn/scikit-learn/discussions
គេហទំព័រ៖ https://scikit-learn.org
LinkedIn៖ https://www.linkedin.com/company/scikit-learn
យូធូប៖ https://www.youtube.com/channel/UCJosFjYm0ZYVUARxuOZqnnw/playlists
ហ្វេសប៊ុក៖ https://www.facebook.com/scikitlearnofficial/
Instagram៖ https://www.instagram.com/scikitlearnofficial/
TikTok៖ https://www.tiktok.com/@scikit.learn
ការដកស្រង់
ប្រសិនបើអ្នកប្រើ scikit-learn នៅក្នុងការបោះពុម្ភផ្សាយបែបវិទ្យាសាស្ត្រ យើងនឹងពេញចិត្តចំពោះការដកស្រង់៖ https://scikit-learn.org/stable/about.html#citing-scikit-learn